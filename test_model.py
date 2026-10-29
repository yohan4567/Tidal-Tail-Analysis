import json
import unittest
from pathlib import Path
import numpy as np
from stream_model import validate, mass_weights
from fit_rate import fit


class ScientificChecks(unittest.TestCase):
    def test_total_mass_both_arms(self):
        for n in (20,2000):
            w=mass_weights(n,1000,2)
            self.assertAlmostEqual(w.sum(),2000)
            self.assertAlmostEqual(w[:n//2].sum(),1000)

    def test_expected_poisson_fit(self):
        s=np.array([1,3,8,3,1.]); b=np.ones(5)*4
        r=fit(3*s+2*b,s,b)
        self.assertAlmostEqual(r['amplitude'],3,places=4)
        self.assertAlmostEqual(r['background_scale'],2,places=4)

    def test_degenerate_and_invalid_data(self):
        with self.assertRaises(ValueError): fit([2,2,2],[1,1,1],[2,2,2])
        with self.assertRaises(ValueError): fit([1.2,2,3],[1,2,3],[3,2,1])
        with self.assertRaises(ValueError): fit([[1,2,3]],[1,2,3],[3,2,1])
        with self.assertRaises(ValueError): fit([100000.1,2,3],[1,2,3],[3,2,1])

    def test_config(self):
        c=json.loads((Path(__file__).parent/'configs/Pal5.json').read_text())
        validate(c)
        c['mass_msun']=-1
        with self.assertRaises(ValueError): validate(c)


if __name__=='__main__': unittest.main()
