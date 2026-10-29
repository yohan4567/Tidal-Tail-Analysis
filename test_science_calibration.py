import unittest
from unittest.mock import patch
import io
import tempfile
from pathlib import Path
import numpy as np
import pandas as pd

from population_model import build_population
from prepare_science_calibration import validate_quality,spatial_parent
from selection_bridge import DEFAULT_STUDY, load_study
from selection_model import calibrate_quality,download_parent
from test_observation import parent
from fetch_isochrone import Defaults
from fetch_parent_tiles import tile_query
from response_model import ErrorLibrary


class ScienceCalibrationChecks(unittest.TestCase):
    def test_escape_tilt_changes_population_but_preserves_mass(self):
        iso=pd.DataFrame(dict(Mini=[.1,.5,.8],Mass=[.1,.5,.7],int_IMF=[0.,1.,1.2],Gmag=[12,6,3],G_BPmag=[13,7,4],G_RPmag=[11,5,2]))
        def make(tilt):
            return build_population(iso,mass_basis='isochrone_current_stellar_mass',source_description='test',age_years=1e10,metallicity_dex=-1,synthetic=True,escape_mass_tilt=tilt)
        a,b=make(-1),make(1)
        for pop in (a,b): self.assertAlmostEqual(np.sum(pop.table.current_mass_msun*pop.table.number_per_msun),1.)
        self.assertGreater(a.table.number_per_msun.iloc[0]/a.table.number_per_msun.iloc[-1],b.table.number_per_msun.iloc[0]/b.table.number_per_msun.iloc[-1])

    def test_real_parent_cannot_hide_a_pm_selection(self):
        T=load_study(DEFAULT_STUDY,True)
        with self.assertRaisesRegex(ValueError,'member/quality'):
            calibrate_quality(parent(),T,center=(10,0),axes=[[-1,1],[-1,1],[17,19],[0,1]],
                parent_metadata=dict(parent_is_unfiltered=True,source_description='test',query_adql='SELECT * FROM gaiadr3.gaia_source WHERE pmra > 0'))

    def test_quality_holdout_counts_missing_sources(self):
        T=load_study(DEFAULT_STUDY,True); train=parent(100); test=parent(100)
        train.loc[:19,'bp_rp']=np.nan; test.loc[:19,'bp_rp']=np.nan
        meta=dict(parent_is_unfiltered=True,source_description='test',synthetic_calibration=True)
        q=calibrate_quality(train,T,center=(10,0),axes=[[-1,1],[-1,1],[17,19],[0,1]],parent_metadata=meta)
        stats,cells=validate_quality(train,test,q,T,(10,0))
        self.assertEqual(stats['available']['n_holdout'],100)
        self.assertAlmostEqual(stats['available']['observed_success_fraction'],.8)
        self.assertEqual(stats['quality']['n_holdout'],80)
        self.assertEqual(len(cells),2)

    def test_cmd_defaults_keep_only_selected_options(self):
        form=Defaults(); form.feed('<input name="x" value="a" type="radio"><input name="x" value="b" type="radio" checked><select name="s"><option value="c"><option value="d" selected></select>')
        self.assertEqual(form.values,dict(x='b',s='d'))

    def test_dark_termination_is_not_a_remnant_census(self):
        iso=pd.DataFrame(dict(Mini=[.1,.8,1.2],Mass=[.1,.7,.5],int_IMF=[0.,1.,2.],Gmag=[12,3,30],G_BPmag=[13,4,31],G_RPmag=[11,2,29],label=[1,3,9]))
        pop=build_population(iso,mass_basis='isochrone_current_stellar_mass',source_description='test',age_years=1e10,metallicity_dex=-1,synthetic=True)
        self.assertEqual(len(pop.table),2)
        self.assertEqual(pop.metadata['excluded_dark_termination_nodes'],1)
        self.assertAlmostEqual(np.sum(pop.table.current_mass_msun*pop.table.number_per_msun),1.)

    def test_indexed_random_parent_has_no_quality_cut(self):
        q=tile_query([10],6,12,21,10)
        where=q.split('WHERE')[1]
        fac=2**35*4**6
        self.assertIn(f'source_id BETWEEN {10*fac} AND {11*fac-1}',where)
        self.assertIn('MOD(random_index,10)',where)
        T=load_study(DEFAULT_STUDY,True)
        meta=dict(parent_is_unfiltered=True,source_description='test',query_adql=q)
        calibrate_quality(parent(),T,center=(10,0),axes=[[-1,1],[-1,1],[17,19],[0,1]],parent_metadata=meta)
        meta['query_adql']=q+' AND phot_g_mean_flux_over_error > 10'
        with self.assertRaisesRegex(ValueError,'unsupported cuts'):
            calibrate_quality(parent(),T,center=(10,0),axes=[[-1,1],[-1,1],[17,19],[0,1]],parent_metadata=meta)

    def test_interrupted_parent_is_never_exported(self):
        class Broken(io.BytesIO):
            def read1(self,n=-1): raise EOFError('interrupted transfer')
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'parent.csv.gz'
            with patch('selection_model.open_url',return_value=Broken()):
                with self.assertRaises(EOFError): download_parent('test query',out)
            self.assertFalse(out.exists())
            self.assertFalse(out.with_suffix('.json').exists())
            self.assertFalse(out.with_suffix(out.suffix+'.download.vot').exists())

    def test_prepared_donor_support_preserves_seeded_joint_errors(self):
        T=load_study(DEFAULT_STUDY,True); d=parent(100)
        d['pmra_error']=np.linspace(.01,.2,len(d))
        errors=ErrorLibrary(d,T,dict(parent_is_unfiltered=True,source_description='test'),{},(10,0))
        g=np.full(5,18.); c=np.full(5,.5); x=y=np.zeros(5)
        support=errors.donor_support(g,c,x,y)
        direct=errors.sample(g,c,x,y,np.random.default_rng(42))
        prepared=errors.sample(g,c,x,y,np.random.default_rng(42),support=support)
        for a,b in zip(direct,prepared): np.testing.assert_array_equal(a,b)

    def test_aip_field_id_does_not_replace_source_id_name(self):
        data=b'''<?xml version="1.0"?><VOTABLE version="1.3" xmlns="http://www.ivoa.net/xml/VOTable/v1.3"><RESOURCE type="results"><INFO name="QUERY_STATUS" value="OK"/><TABLE><FIELD ID="datalinkid" name="source_id" datatype="long"/><DATA><TABLEDATA><TR><TD>123456789012345678</TD></TR></TABLEDATA></DATA></TABLE></RESOURCE></VOTABLE>'''
        with tempfile.TemporaryDirectory() as tmp:
            out=Path(tmp)/'parent.csv.gz'
            with patch('selection_model.open_url',return_value=io.BytesIO(data)):
                download_parent('test query',out)
            d=pd.read_csv(out)
            self.assertIn('source_id',d)
            self.assertEqual(int(d.source_id.iloc[0]),123456789012345678)

    def test_parent_sky_masks_exclude_crowded_cores_only(self):
        T=load_study(DEFAULT_STUDY,True); d=parent(4); d.ra=[10.,10.5,12.2,13.]
        report=dict(center_ra_deg=10.,center_dec_deg=0.,tidal_radius_deg=.1,neighbour_masks=[['other',2.,0.,.3]])
        out=spatial_parent(d,T,report)
        self.assertEqual(out.index.tolist(),[1,3])


if __name__=='__main__': unittest.main()
