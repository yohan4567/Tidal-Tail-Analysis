"""Copy verified completed simulations or checkpoints into a different parameter grid."""
import hashlib,json,shutil
from pathlib import Path
from importlib.metadata import version
import numpy as np
from provenance import sha256_file
HERE=Path(__file__).resolve().parent


def reuse_simulation(source,destination,config):
    source,destination=map(Path,(source,destination))
    if destination.exists():raise ValueError('Reuse destination already exists')
    signature_path=source/'simulation_input.json'
    if not signature_path.exists():return None
    signature=json.loads(signature_path.read_text())
    if signature['config']!=config:raise ValueError('Reuse physical/numerical configuration mismatch')
    for name,digest in signature['source_sha256'].items():
        if sha256_file(HERE/name)!=digest:raise ValueError('Reuse numerical source mismatch: '+name)
    for name,value in signature['versions'].items():
        if version(name)!=value:raise ValueError('Reuse package version mismatch: '+name)
    complete=(source/'run.json').exists();checkpoint=source/'integration_checkpoint.npz';count=0
    if complete:
        metadata=json.loads((source/'run.json').read_text())
        if metadata['config']!=config or sha256_file(source/'particles.csv')!=metadata['particles_sha256']:
            raise ValueError('Reuse completed output checksum/config mismatch')
        count=config['n_particles']
    elif checkpoint.exists():
        model=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
        with np.load(checkpoint,allow_pickle=False) as data:
            if json.loads(str(data['identity'].item()))['model']!=model:raise ValueError('Reuse checkpoint model mismatch')
            states=data['final_cartesian']
            if states.ndim!=2 or states.shape[1]!=6 or len(states)>config['n_particles'] or not np.isfinite(states).all():
                raise ValueError('Reuse invalid checkpoint')
            count=len(states)
    files={p.name:sha256_file(p) for p in source.iterdir() if p.is_file() and not p.name.endswith('.tmp.npz')}
    shutil.copytree(source,destination,ignore=shutil.ignore_patterns('*.tmp.npz'))
    for name,digest in files.items():
        if sha256_file(source/name)!=digest or sha256_file(destination/name)!=digest:
            raise ValueError('Source changed during reuse or copied checksum mismatch')
    audit=dict(source=str(source.resolve()),destination=str(destination.resolve()),complete=complete,
               retained_tracers=count,files_sha256=files,validation='Exact configuration, numerical sources, package versions and file checksums verified')
    (destination.parent/'reuse_provenance.json').write_text(json.dumps(audit,indent=2))
    return audit
