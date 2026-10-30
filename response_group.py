"""Compute compatible population/quality responses together, with native exports."""
import argparse
import json
from pathlib import Path
import numpy as np
from response_model import run, DEFAULT_STUDY
from population_model import Population
from response_shared import expected_shared
from validate_science import check_response
from provenance import sha256_file


class GroupProvider:
    def __init__(self,populations):
        self.populations=populations;self.results=None;self.diagnostics=None

    def __call__(self,particles,cases,*args,diagnostics,**options):
        index=next((i for i,p in enumerate(self.populations) if p.table.equals(cases[0]['population'].table)),None)
        if index is None:raise ValueError('Population not registered in response group')
        count=len(cases)
        if self.results is None:
            all_cases=[dict(population=p,quality=c['quality']) for p in self.populations for c in cases]
            self.diagnostics=[dict(fitmask=diagnostics[0]['fitmask']) for _ in all_cases]
            self.results=expected_shared(particles,all_cases,*args,diagnostics=self.diagnostics,**options)
        lo=index*count;hi=lo+count
        for target,source in zip(diagnostics,self.diagnostics[lo:hi]):target.update(source)
        return self.results[lo:hi]


def compatible(specs,populations):
    # All sampling and observation settings must agree; descriptions and population
    # weights may differ. Paths are resolved before this comparison.
    ignored={'population','scenario','source_description','cache_dir'}
    clean=lambda s:{k:v for k,v in s.items() if k not in ignored}
    if any(clean(s)!=clean(specs[0]) for s in specs[1:]):return False
    table=populations[0].table.drop(columns='number_per_msun')
    support=populations[0].table.number_per_msun.to_numpy()>0
    return all(p.metadata['mass_basis']==populations[0].metadata['mass_basis'] and
        p.table.drop(columns='number_per_msun').equals(table) and
        np.array_equal(p.table.number_per_msun.to_numpy()>0,support) for p in populations)


def run_group(selection,simulation,configs,outputs,study=DEFAULT_STUDY):
    configs,outputs=list(map(Path,configs)),list(map(Path,outputs))
    if not configs or len(configs)!=len(outputs):raise ValueError('Matching nonempty config/output lists required')
    if len({p.resolve() for p in outputs})!=len(outputs):raise ValueError('Duplicate response destinations')
    pending=[];results=[None]*len(configs)
    for i,(cfg,out) in enumerate(zip(configs,outputs)):
        if (out/'response.json').exists():
            check_response(out)
            meta=json.loads((out/'response.json').read_text(encoding='utf-8'))
            if meta['input_sha256'].get(str(cfg.resolve()))!=sha256_file(cfg):raise ValueError('Stale response config')
            results[i]=meta
        else:
            if out.exists():raise ValueError('Incomplete output directory must be inspected before retry: '+str(out))
            pending.append(i)
    specs=[];populations=[]
    for i in pending:
        cfg=configs[i];s=json.loads(cfg.read_text(encoding='utf-8'));base=cfg.resolve().parent
        for k in ('population','quality','cache_dir','rate_history'):
            if s.get(k):s[k]=str((base/s[k]).resolve())
        for k in ('path','metadata'):s['errors'][k]=str((base/s['errors'][k]).resolve())
        for k in ('path','metadata'):
            if s['completeness'].get(k):s['completeness'][k]=str((base/s['completeness'][k]).resolve())
        specs.append(s);populations.append(Population.load(s['population']))
    provider=GroupProvider(populations) if pending and compatible(specs,populations) else None
    for i in pending:
        results[i]=run(selection,simulation,configs[i],outputs[i],study,backend='shared',_provider=provider)
        check_response(outputs[i])
        print('Response export complete:',outputs[i],flush=True)
    return results


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--selection',type=Path,required=True);p.add_argument('--run',type=Path,required=True)
    p.add_argument('--configs',type=Path,nargs='+',required=True);p.add_argument('--outputs',type=Path,nargs='+',required=True)
    p.add_argument('--study',type=Path,default=DEFAULT_STUDY);a=p.parse_args()
    run_group(a.selection,a.run,a.configs,a.outputs,a.study)
