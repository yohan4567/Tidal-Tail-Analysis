import unittest
from types import SimpleNamespace
import numpy as np
import pandas as pd
from population_model import Population
from response_model import expected_template
from response_shared import expected_shared

class Cuts:
    G_MAX=20.5
    def gnomonic(self,ra,dec,ra0,dec0):return np.asarray(ra)-ra0,np.asarray(dec)-dec0
    def plx_select(self,d,distance):return abs(d['parallax']-1/distance)<3*d['parallax_error']+.05
    def pm_select(self,d,pm,distance,sig0):return ((d['pmra']-pm[0])**2+(d['pmdec']-pm[1])**2)<.15
    def cmd_select(self,d,mids,med,wid,dg):
        g=d['phot_g_mean_mag']-(0 if dg is None else dg)
        return (g>=mids.min()-.3)&(g<=mids.max()+.3)&(abs(d['bp_rp']-np.interp(g,mids,med))<np.interp(g,mids,wid))
    def orbit_expectation(self,row,R,x,y):return (row.pmra+.01*x,row.pmdec+.01*y),row.R_sun+.02*x,0.,0.

class Errors:
    max_g_error=.08
    def __init__(self):self.calls=0
    def donor_support(self,g,c,x,y):return None
    def sample(self,g,c,x,y,rng,support=None):
        self.calls+=1;rng.random(len(g))
        root=np.array([[.02,0,0],[.01,.03,0],[0,0,.03]])
        return np.broadcast_to(root,(len(g),3,3)),np.full((len(g),3),.08)

class Survey:
    def query(self,ra,dec,g,x,y):return np.where(x>.5,0.,.9)

class Quality:
    def __init__(self,p):self.p=p
    def query(self,x,y,g,c):return np.full(len(g),.8),np.full(len(g),self.p)

class SharedResponseTests(unittest.TestCase):
    def fixture(self,drift=False):
        table=pd.DataFrame(dict(M_G=[3.,4.,8.],M_BP=[3.3,4.3,8.3],M_RP=[2.5,3.5,7.5],current_mass_msun=[.5,.7,.9],number_per_msun=[1.,.2,.3]))
        pops=[]
        for f in [1.,.7,1.4]:
            t=table.copy();t.number_per_msun*=np.array([1.,f,1./f]);t.number_per_msun/=np.sum(t.current_mass_msun*t.number_per_msun)
            pops.append(Population(t,dict(mass_basis='isochrone_current_stellar_mass',source_description='Software fixture')))
        cases=[dict(population=pop,quality=Quality(q)) for pop in pops for q in [.2,.9,0.]]
        stars=pd.DataFrame(dict(ra_deg=[9.7,10.,10.2,10.7,13.],dec_deg=[0.,.2,0.,0.,0.],distance_kpc=[10.]*5,pmra_masyr=[1.]*5,pmdec_masyr=[2.]*5))
        row=SimpleNamespace(ra=10.,dec=0.,R_sun=10.,pmra=1.,pmdec=2.)
        P=dict(R=1.,sig0=5.,drift=dict(corrected=drift),mids=np.array([17.,18.,19.]),med=np.full(3,.8),wid=np.full(3,.3),masks=[])
        meta=dict(radius_deg=1.,pm0=[1.,2.],pm_box=2.,g_limit=20.5)
        return stars,cases,Survey(),row,P,meta
    def test_all_outputs_identical_and_noise_draws_shared(self):
        for drift in [False,True]:
            stars,cases,survey,row,P,meta=self.fixture(drift);T=Cuts();inside=np.ones((10,10),bool);edges=np.linspace(-1,1,11)
            opts=dict(a_g=0.,a_bp=0.,a_rp=0.,weights=np.array([2.,1.,0.,1.,1.]),draws=8,seed=123,max_pairs=6)
            refs=[];refdiag=[];old=Errors()
            for case in cases:
                d=dict(fitmask=inside);refdiag.append(d)
                refs.append(expected_template(stars,case['population'],case['quality'],survey,old,T,row,P,meta,inside,edges,diagnostics=d,**opts))
            new=Errors();diags=[dict(fitmask=inside) for _ in cases]
            shared=expected_shared(stars,cases,survey,new,T,row,P,meta,inside,edges,diagnostics=diags,**opts)
            for ref,value,rd,d in zip(refs,shared,refdiag,diags):
                np.testing.assert_array_equal(ref[0],value[0]);np.testing.assert_array_equal(ref[1],value[1]);self.assertEqual(ref[2],value[2])
                np.testing.assert_array_equal(rd['tracer_selected_counts'],d['tracer_selected_counts'])
            self.assertEqual(old.calls,new.calls*len(cases))
    def test_rejects_different_nodes_and_zero_weight_support(self):
        stars,cases,survey,row,P,meta=self.fixture();opts=dict(a_g=0,a_bp=0,a_rp=0,weights=np.ones(5),draws=4)
        cases[-1]['population'].table.loc[0,'M_G']+=.1
        with self.assertRaises(ValueError):expected_shared(stars,cases,survey,Errors(),Cuts(),row,P,meta,np.ones((2,2),bool),np.linspace(-1,1,3),**opts)
        stars,cases,survey,row,P,meta=self.fixture();cases[-1]['population'].table.loc[0,'number_per_msun']=0
        with self.assertRaises(ValueError):expected_shared(stars,cases,survey,Errors(),Cuts(),row,P,meta,np.ones((2,2),bool),np.linspace(-1,1,3),**opts)

if __name__=='__main__':unittest.main()
