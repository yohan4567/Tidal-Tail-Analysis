import copy,unittest
from pathlib import Path
from recover_compact_repairs import rebase,semantic_spec,HERE

class CompactRepairTests(unittest.TestCase):
    def test_only_whole_path_prefix_is_relocated(self):
        old=Path('old').resolve();new=Path('short').resolve()
        value={'input':str(old/'data/file.json'),'description':'preserve original science',
               'other':str(old)+'_another/file.json','number':.5}
        moved=rebase(value,old,new)
        self.assertEqual(moved['input'],str(new/'data/file.json'))
        self.assertEqual(moved['description'],value['description'])
        self.assertEqual(moved['other'],value['other']);self.assertEqual(moved['number'],.5)
    def test_relocation_cannot_relax_donor_distance_or_sampling(self):
        spec={'population':'old','quality':'old','cache_dir':'old',
              'errors':{'path':'old','metadata':'old','max_feature_distance':2.,'neighbours':20},
              'observation_draws':32,'calibration_draws':8,'seed':42}
        moved=copy.deepcopy(spec);moved.update(population='new',quality='new',cache_dir='new')
        moved['errors'].update(path='new',metadata='new')
        self.assertEqual(semantic_spec(spec),semantic_spec(moved))
        moved['errors']['max_feature_distance']=4.
        self.assertNotEqual(semantic_spec(spec),semantic_spec(moved))
    def test_compact_cache_fits_conservative_windows_budget(self):
        path=HERE/'outputs/rc2/NGC6341'/('f'*64)/'age_contributions.csv'
        self.assertLess(len(str(path)),248)

if __name__=='__main__':unittest.main()
