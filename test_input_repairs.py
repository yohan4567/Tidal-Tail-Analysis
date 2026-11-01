"""Scientific measure preservation, exact parent coverage, and safe I/O retry checks."""
import json,re,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
import pandas as pd
from population_quadrature_repair import mass_quadrature
from checkpoint_io import replace_retry,retry_path_replaces
from adaptive_parent import fetch,split_patch
from provenance import sha256_file


class QuadratureTests(unittest.TestCase):
    def test_unique_measure_matches_analytic_integral(self):
        m=np.array([.1,.2,.4,.8])
        for slope in (-2.35,-1.,0.,1.):
            w,a=mass_quadrature(m,slope,np.ones(4))
            target=np.log(8) if slope==-1 else (.8**(slope+1)-.1**(slope+1))/(slope+1)
            self.assertAlmostEqual(w.sum(),target)
            self.assertEqual(a['tied_nodes'],0)

    def test_tied_nodes_retain_measure_and_supplied_phase_ratio(self):
        m=np.array([.1,.1,.1,.2,.4,.4]);prior=np.array([1,2,3,1,4,2.])
        w,a=mass_quadrature(m,-1.,prior)
        unique,_=mass_quadrature(np.unique(m),-1.,np.ones(3))
        np.testing.assert_allclose([w[:3].sum(),w[3],w[4:].sum()],unique)
        np.testing.assert_allclose(w[:3]/w[0],[1,2,3])
        self.assertAlmostEqual(w[4]/w[5],2.)
        self.assertTrue((w>0).all());self.assertEqual(a['tied_nodes'],3)
        current=np.array([.1,.1,.1,.2,.3,.25])
        normalized=w/np.sum(w*current)
        self.assertAlmostEqual(np.sum(normalized*current),1.)

    def test_unsupported_coordinates_rejected(self):
        for masses in ([.1,.1],[.2,.1],[0.,.2],[np.nan,.2]):
            with self.assertRaises(ValueError):mass_quadrature(masses,-1,np.ones(2))


class CheckpointTests(unittest.TestCase):
    def test_retry_then_real_atomic_replace(self):
        with tempfile.TemporaryDirectory() as d:
            source=Path(d)/'new';dest=Path(d)/'old';source.write_text('new');dest.write_text('old')
            calls=[]
            def operation(p,q):
                calls.append(1)
                if len(calls)<3:raise PermissionError('temporary lock')
                return p.replace(q)
            replace_retry(source,dest,operation=operation,sleep=lambda _:None)
            self.assertEqual(dest.read_text(),'new');self.assertEqual(len(calls),3)

    def test_permanent_lock_is_bounded(self):
        calls=[]
        def locked(*args):calls.append(1);raise PermissionError('locked')
        with self.assertRaises(PermissionError):replace_retry('a','b',attempts=3,operation=locked,sleep=lambda _:None)
        self.assertEqual(len(calls),3)
        with self.assertRaises(FileNotFoundError):replace_retry('missing-file','b')

    def test_context_restores_even_on_error(self):
        original=Path.replace
        with self.assertRaises(RuntimeError):
            with retry_path_replaces():
                self.assertIsNot(Path.replace,original);raise RuntimeError('test')
        self.assertIs(Path.replace,original)

    def test_real_integrator_recovers_from_locked_checkpoint(self):
        from test_batch_dynamics import BatchDynamicsTests
        from batch_dynamics import integrate_released_batched
        args=BatchDynamicsTests().fixture();original=Path.replace;calls=[]
        def locked_once(source,destination):
            calls.append(1)
            if len(calls)==1:raise PermissionError('injected transient file lock')
            return original(source,destination)
        with tempfile.TemporaryDirectory() as d:
            checkpoint=Path(d)/'state.npz'
            with patch.object(Path,'replace',locked_once),retry_path_replaces():
                result=integrate_released_batched(*args,batch_size=2,checkpoint_path=checkpoint,model_identity='repair_test',checkpoint_every=1)
            with retry_path_replaces():
                resumed=integrate_released_batched(*args,batch_size=2,checkpoint_path=checkpoint,model_identity='repair_test',checkpoint_every=1)
            np.testing.assert_array_equal(resumed,result)
            self.assertGreater(len(calls),1)


class ParentTests(unittest.TestCase):
    def test_overflow_recursion_has_exact_nonduplicated_coverage_and_resume(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);calls=[]
            def download(query,out,maxrec,endpoint):
                ranges=re.findall(r'source_id BETWEEN (\d+) AND (\d+)',query)
                widths=[int(b)-int(a)+1 for a,b in ranges];calls.append(query)
                # Require both splitting multi-tile batches and refining a tile.
                if len(ranges)>1 or widths[0]>2**35*4**5:raise RuntimeError('TAP incomplete (OVERFLOW)')
                rows=[dict(source_id=int(a)+10,ra=20.,dec=0.) for a,b in ranges]
                pd.DataFrame(rows).to_csv(out,index=False,compression='infer')
                meta=dict(query_status='OK',query_adql=query,sha256=sha256_file(out),endpoint=endpoint)
                out.with_suffix('.json').write_text(json.dumps(meta))
            out=root/'first.csv.gz'
            result=fetch(out,20.,0.,1.,pixels=[123,124],downloader=download)
            data=pd.read_csv(out);self.assertEqual(len(data),8);self.assertFalse(data.source_id.duplicated().any())
            self.assertEqual(len(result['leaf_patches']),8)
            count=len(calls)
            second=fetch(root/'second.csv.gz',20.,0.,1.,pixels=[123,124],downloader=download,reuse_dirs=[root/'first.csv.gz.parts'])
            self.assertTrue(all(p['reused'] for p in second['leaf_patches']))
            # Parent overflow requests recur, but successful leaves are never downloaded twice.
            self.assertEqual(len(calls)-count,3)
            pd.testing.assert_frame_equal(data,pd.read_csv(root/'second.csv.gz'))

    def test_checksum_corruption_blocks_cache(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);cached=root/'cache';cached.mkdir()
            from fetch_parent_tiles import tile_query
            file=cached/'one.csv.gz';file.write_bytes(b'corrupt')
            file.with_suffix('.json').write_text(json.dumps(dict(query_status='OK',query_adql=tile_query([123],6,12.,21.,10),sha256='bad')))
            with self.assertRaisesRegex(ValueError,'checksum'):
                fetch(root/'out.csv.gz',20.,0.,1.,pixels=[123],reuse_dirs=[cached])
            self.assertFalse((root/'out.csv.gz').exists())

    def test_finest_overflow_rejected(self):
        with self.assertRaises(RuntimeError):split_patch([3],12)


class ProgressTests(unittest.TestCase):
    def test_external_pc_exclusion_survives_preparation_status_updates(self):
        from dispatch_repair_queue import selected_clusters
        available=('NGC3201','NGC362','NGC6121');overrides={'excluded_clusters':['NGC362']}
        for rewritten_stage in ('downloading_or_calibrating','ready_for_conditional_inference'):
            # Dispatch consults separate overrides, never the producer's stage.
            self.assertEqual(selected_clusters(available,overrides),('NGC3201','NGC6121'))
        with self.assertRaises(ValueError):selected_clusters(available,{'excluded_clusters':['unknown']})

    def test_repair_status_is_shown_separately_from_original_pass(self):
        from refresh_cluster_progress import refresh
        with tempfile.TemporaryDirectory() as d:
            base=Path(d);root=base/'original';root.mkdir();repair=base/'input_repairs_test';repair.mkdir()
            (repair/'repair_manifest.json').write_text(json.dumps(dict(source=str(root))))
            (repair/'input_status.json').write_text(json.dumps(dict(NGC362=dict(stage='ready_for_conditional_inference'))))
            # Dry-run reporting: no atomic replacements or live progress writes.
            with patch.object(Path,'write_text'),patch.object(Path,'replace'):
                report=refresh(root)
            entry=next(e for e in report['entries'] if e['cluster']=='NGC362')
            self.assertEqual(entry['input_repair']['stage'],'ready_for_conditional_inference')
            self.assertEqual(entry['status'],'queued')
            self.assertNotIn('input_repair',next(e for e in report['entries'] if e['cluster']=='Pal5'))


if __name__=='__main__':unittest.main()
