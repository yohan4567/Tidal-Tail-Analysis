"""Explicit spatial pooling of sparse quality cells, with shared beta draws.

Availability retains its original sky/magnitude grid. Quality retains magnitude
and colour; a sparse magnitude/colour stratum pools all its spatial cells.
This is a conditional exchangeability assumption, not added calibration data.
"""
from copy import deepcopy
import numpy as np
from selection_model import QualityMap, lookup_indices


class SpatialPooledQualityMap(QualityMap):
    def __init__(self, arrays, metadata, draw_seed=None):
        super().__init__(arrays, metadata, draw_seed)
        if metadata.get('quality_pooling', {}).get('policy') != 'spatial_stratum_if_sparse_v1':
            raise ValueError('Explicit spatial pooling metadata required')
        ids=np.asarray(arrays['quality_group_ids'])
        n=np.asarray(arrays['quality_group_n']);k=np.asarray(arrays['quality_group_k'])
        if ids.shape!=arrays['n_available'].shape or ids.dtype.kind not in 'iu' or ids.min()<0:
            raise ValueError('Invalid quality group IDs')
        if n.ndim!=1 or k.shape!=n.shape or ids.max()>=len(n) or set(np.unique(ids))!=set(range(len(n))):
            raise ValueError('Invalid quality group counts')
        for group in range(len(n)):
            mask=ids==group
            if n[group]!=arrays['n_available'][mask].sum() or k[group]!=arrays['k_quality'][mask].sum():
                raise ValueError('Pooled counts disagree with original cells')
        if (n<self.minimum).any():
            raise ValueError('Spatial pooling still unsupported; obtain more parent data')
        self.group_ids=ids;self.group_n=n
        if draw_seed is None:
            self.group_probability=(k+1)/(n+2)
        else:
            rng=np.random.default_rng(draw_seed)
            # Retain the native availability draw and use the subsequent stream
            # once per disjoint quality group, shared by every member cell.
            self.probabilities['available']=rng.beta(arrays['k_available']+1,
                arrays['n_parent']-arrays['k_available']+1)
            self.group_probability=rng.beta(k+1,n-k+1)
        self.probabilities['quality']=self.group_probability[ids]

    def query(self,x,y,g,g_rp):
        i3,ok3=lookup_indices(x,y,g,self.axes[:3])
        i4,ok4=lookup_indices(x,y,g,self.axes,color=g_rp)
        if not (ok3&ok4).all():
            raise ValueError('Quality query outside calibrated spatial/magnitude/colour grid')
        if (self.arrays['n_parent'][i3]<self.minimum).any() or (self.group_n[self.group_ids[i4]]<self.minimum).any():
            raise ValueError('Quality calibration unsupported after spatial pooling')
        return self.probabilities['available'][i3],self.probabilities['quality'][i4]


def pool_sparse_spatial_cells(quality,support_counts=None):
    """Partition using eligible counts only, never successes or a fitted rate."""
    arrays={k:v.copy() for k,v in quality.arrays.items() if not k.startswith('quality_group_')}
    n=arrays['n_available'];k=arrays['k_quality'];minimum=quality.minimum
    support=n if support_counts is None else np.asarray(support_counts)
    if support.shape!=n.shape or (support<0).any() or (support>n).any():
        raise ValueError('Invalid independent-training support counts')
    ids=np.full(n.shape,-1,dtype=np.int64);group_n=[];group_k=[];pooled=[]
    for ig in range(n.shape[2]):
        for ic in range(n.shape[3]):
            if (support[:,:,ig,ic]<minimum).any():
                index=len(group_n);ids[:,:,ig,ic]=index
                group_n.append(int(n[:,:,ig,ic].sum()));group_k.append(int(k[:,:,ig,ic].sum()))
                pooled.append(dict(group=index,g_bin=ig,colour_bin=ic,
                    g_range=arrays['g_edges'][ig:ig+2].tolist(),
                    colour_range=arrays['c_edges'][ic:ic+2].tolist(),n=group_n[-1]))
            else:
                for ix in range(n.shape[0]):
                    for iy in range(n.shape[1]):
                        ids[ix,iy,ig,ic]=len(group_n)
                        group_n.append(int(n[ix,iy,ig,ic]));group_k.append(int(k[ix,iy,ig,ic]))
    arrays.update(quality_group_ids=ids,quality_group_n=np.asarray(group_n,dtype=np.int64),
        quality_group_k=np.asarray(group_k,dtype=np.int64))
    meta=deepcopy(quality.metadata);meta.pop('table_sha256',None)
    meta['quality_pooling']=dict(policy='spatial_stratum_if_sparse_v1',
        selection_rule='Pool all sky cells in a G/colour stratum if any eligible training count is below min_parent; no quality outcomes or rate fits select groups.',
        support_count_basis='all_parent' if support_counts is None else 'independent_training_partition',
        pooled_strata=pooled,groups=len(group_n),original_cells=int(n.size),
        original_low_count_cells=int((n<minimum).sum()),
        assumption='Quality probability is spatially exchangeable within each pooled G/colour stratum; availability retains the original sky/G grid.',
        uncertainty='One Beta(k+1,n-k+1) draw per disjoint group, shared by all its cells and particles.',
        physical_rate_measurement_validated=False)
    return SpatialPooledQualityMap(arrays,meta)


def validate_pooling(full_quality,holdout_cells,screening='gaussian_v1'):
    """Validate quality predictions in the original fine sky/G/colour cells.

    Group definitions use eligible counts. Quality-success counts in the holdout
    are subtracted before fitting the beta predictive distribution.
    """
    if screening not in ('gaussian_v1','exact_beta_binomial_v2'):
        raise ValueError('Unknown held-out screening method')
    from scipy.stats import betabinom,norm
    train_n=full_quality.arrays['n_available'].copy();train_k=full_quality.arrays['k_quality'].copy()
    held=holdout_cells.loc[holdout_cells.stage=='quality']
    for row in held.itertuples():
        index=tuple(map(int,row.cell.split(',')))
        train_n[index]-=int(row.holdout_n);train_k[index]-=int(row.holdout_k)
    if (train_n<0).any() or (train_k<0).any() or (train_k>train_n).any():
        raise ValueError('Held-out counts are inconsistent with the full calibration')
    candidate=pool_sparse_spatial_cells(full_quality,support_counts=train_n)
    rows=[]
    for row in held.itertuples():
        index=tuple(map(int,row.cell.split(',')));group=int(candidate.group_ids[index]);mask=candidate.group_ids==group
        n=int(train_n[mask].sum());k=int(train_k[mask].sum())
        if n<full_quality.minimum:raise ValueError('Insufficient training support in pooled holdout check')
        prob=(k+1)/(n+2);count=int(row.holdout_n)
        variance=count*prob*(1-prob)*(count+n+2)/(n+3)
        z=(row.holdout_k-count*prob)/np.sqrt(variance) if variance>0 else 0.
        tail=min(1.,2*min(betabinom.cdf(int(row.holdout_k),count,k+1,n-k+1),
                          betabinom.sf(int(row.holdout_k)-1,count,k+1,n-k+1)))
        rows.append(dict(cell=row.cell,group=group,training_n=n,holdout_n=count,holdout_k=int(row.holdout_k),
            probability=float(prob),predictive_z=float(z),exact_predictive_two_sided_tail=float(tail)))
    largest=max((abs(row['predictive_z']) for row in rows),default=float('inf'))
    cutoff=float(2*norm.sf(3.));minimum_tail=min((row['exact_predictive_two_sided_tail'] for row in rows),default=0.)
    passed=bool(rows) and (largest<=3 if screening=='gaussian_v1' else minimum_tail>=cutoff)
    return candidate,dict(original_fine_cells_tested=len(rows),max_abs_predictive_z=largest,
        cells_over_abs_z3=sum(abs(row['predictive_z'])>3 for row in rows),
        screening_method=screening,screening_threshold=3.,screening_passed=passed,
        gaussian_screening_passed=bool(rows) and largest<=3,
        min_exact_predictive_two_sided_tail=minimum_tail,exact_tail_threshold=cutoff,
        exact_screening_passed=bool(rows) and minimum_tail>=cutoff,
        exact_tail_definition='min(1, 2*min(P(K<=k),P(K>=k))) under Beta-binomial(n_holdout,k_train+1,n_train-k_train+1); per-cell diagnostic, not a familywise test.',
        interpretation='Held-out field-star check under a stated spatial-pooling assumption; not validation for the escape population or an absolute mass-loss rate.',cells=rows)
