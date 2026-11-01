import json,tempfile,unittest
from io import BytesIO
from pathlib import Path
from urllib.error import HTTPError
from unittest.mock import patch
import repair_parent_support as support
from provenance import sha256_file

class ParentRepairTests(unittest.TestCase):
    selection=dict(center_ra_deg=10.,center_dec_deg=-20.,field_radius_deg=1.)
    def fake_parent(self,query,path,limit,endpoint):
        path=Path(path);path.write_bytes(b'source_id,ra,dec\n')
        m=dict(query_status='OK',parent_is_unfiltered=True,synthetic_calibration=False,
               sha256=sha256_file(path),query_adql=query,n_sources=0,endpoint=endpoint)
        path.with_suffix('.json').write_text(json.dumps(m));return m
    def test_http400_fails_over_without_accepting_error_body(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bright_parent.csv.gz';calls=[]
            def download(q,out,limit,endpoint):
                calls.append(endpoint)
                if endpoint==support.ENDPOINTS[0]:
                    raise HTTPError(endpoint,400,'Bad Request',{},BytesIO(b'invalid TAP query'))
                return self.fake_parent(q,out,limit,endpoint)
            with patch.object(support,'download_parent',side_effect=download),patch.object(support.time,'sleep'):
                m=support.fetch_parent(p,self.selection,bright=True)
            self.assertEqual(calls,[support.ENDPOINTS[0],support.ENDPOINTS[0],support.ENDPOINTS[1]])
            self.assertEqual(m['endpoint'],support.ENDPOINTS[1])
            self.assertEqual(m['random_sampling_modulus'],1)
            self.assertTrue(support.matches(p,10.,-20.,1.05,0.,12.,1))
    def test_clipped_query_is_never_accepted(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'bright.csv.gz'
            with patch.object(support,'download_parent',side_effect=RuntimeError('OVERFLOW')):
                with self.assertRaisesRegex(RuntimeError,'OVERFLOW'):
                    support.fetch_parent(p,self.selection,bright=True)
            self.assertFalse(p.exists())
    def test_server_statement_timeout_switches_to_indexed_tiles(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'parent.csv.gz'
            error=HTTPError('tap',400,'Bad Request',{},BytesIO(b'canceling statement due to statement timeout'))
            def tiled(out,ra,dec,radius,**kwargs):
                m=self.fake_parent('indexed tiles',out,1,'indexed')
                m.update(center=[ra,dec],radius_deg=radius,g_min=12.,g_max=21.,random_sampling_modulus=10)
                Path(out).with_suffix('.json').write_text(json.dumps(m));return m
            with patch.object(support,'download_parent',side_effect=error) as request,patch('adaptive_parent.fetch',side_effect=tiled) as adaptive:
                support.fetch_parent(p,self.selection)
            self.assertEqual(request.call_count,1);adaptive.assert_called_once()
            self.assertEqual(json.loads(p.with_suffix('.json').read_text())['query_adql'],'indexed tiles')
    def test_wrong_geometry_cache_is_not_reused(self):
        with tempfile.TemporaryDirectory() as d:
            old=Path(d)/'old.csv.gz';target=Path(d)/'new.csv.gz'
            self.fake_parent('old',old,1,'old')
            m=json.loads(old.with_suffix('.json').read_text())
            m.update(center=[10.,-20.],radius_deg=2.,g_min=0.,g_max=12.,random_sampling_modulus=1)
            old.with_suffix('.json').write_text(json.dumps(m))
            with patch.object(support,'download_parent',side_effect=self.fake_parent) as request:
                support.fetch_parent(target,self.selection,bright=True,candidates=[old])
            self.assertEqual(request.call_count,1)
    def test_preferred_small_tiles_skip_the_timed_out_cone_strategy(self):
        with tempfile.TemporaryDirectory() as d:
            path=Path(d)/'parent.csv.gz'
            def tiled(out,ra,dec,radius,**kwargs):
                meta=self.fake_parent('small indexed patches',out,1,'indexed')
                meta.update(center=[ra,dec],radius_deg=radius,g_min=12.,g_max=21.,random_sampling_modulus=10)
                Path(out).with_suffix('.json').write_text(json.dumps(meta));return meta
            with patch.object(support,'download_parent') as direct,patch('adaptive_parent.fetch',side_effect=tiled) as adaptive:
                support.fetch_parent(path,self.selection,prefer_tiles=True,tile_level=7,tiles_per_request=1,max_requests=1200)
            direct.assert_not_called()
            self.assertEqual(adaptive.call_args.kwargs['level'],7)
            self.assertEqual(adaptive.call_args.kwargs['tiles_per_request'],1)
            self.assertEqual(adaptive.call_args.kwargs['max_requests'],1200)
    def test_unsupported_checkpoint_is_preserved_and_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'parent.csv.gz';p.write_bytes(b'incomplete')
            with patch.object(support,'download_parent') as request:
                with self.assertRaisesRegex(ValueError,'checkpoint'):support.fetch_parent(p,self.selection)
            self.assertEqual(p.read_bytes(),b'incomplete');request.assert_not_called()
    def test_merge_wait_preserves_verified_tiles_for_restart(self):
        import re
        import pandas as pd
        from adaptive_parent import fetch
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'parent.csv.gz';calls=[]
            def download(query,out,limit,endpoint):
                calls.append(query)
                source_id=int(re.search(r'source_id BETWEEN (\d+)',query).group(1))+1
                pd.DataFrame([dict(source_id=source_id,ra=10.,dec=-20.)]).to_csv(out,index=False,compression='gzip')
                Path(out).with_suffix('.json').write_text(json.dumps(dict(query_status='OK',query_adql=query,
                    sha256=sha256_file(out),endpoint=endpoint)))
            def unavailable():raise RuntimeError('Insufficient RAM for merge')
            with self.assertRaisesRegex(RuntimeError,'RAM'):
                fetch(p,10.,-20.,1.,pixels=[123],downloader=download,before_merge=unavailable)
            self.assertFalse(p.exists())
            result=fetch(p,10.,-20.,1.,pixels=[123],downloader=download,before_merge=lambda:None)
            self.assertEqual(len(calls),1);self.assertTrue(result['leaf_patches'][0]['reused'])

if __name__=='__main__':unittest.main()
