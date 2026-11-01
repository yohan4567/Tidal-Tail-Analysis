"""Bounded retries for transient Windows file locks; numerical code is unchanged."""
from contextlib import contextmanager
from datetime import datetime,timezone
import json,sys,time
from pathlib import Path
from provenance import sha256_file

_replace=Path.replace


def replace_retry(source,destination,*,attempts=12,sleep=time.sleep,operation=None):
    if attempts<1:raise ValueError('Positive retry count required')
    operation=operation or _replace
    for attempt in range(attempts):
        try:return operation(Path(source),destination)
        except PermissionError:
            if attempt==attempts-1:raise
            sleep(min(.15*2**attempt,2.))


@contextmanager
def retry_path_replaces():
    original=Path.replace
    def guarded(path,destination):return replace_retry(path,destination,operation=original)
    Path.replace=guarded
    try:yield
    finally:Path.replace=original


if __name__=='__main__':
    import stream_model
    # Preserve the scientific CLI and its original numerical signatures. The
    # additional audit describes I/O retries; it does not relabel checkpoints.
    args=sys.argv[1:]
    out=Path(args[args.index('--out')+1]).resolve() if '--out' in args else None
    status=dict(policy='bounded_permission_retry_v1',attempts=12,numerical_sources_unchanged=True,
        wrapper_sha256=sha256_file(Path(__file__)),started_utc=datetime.now(timezone.utc).isoformat())
    try:
        with retry_path_replaces():stream_model.main()
        status['stage']='completed'
    except Exception as exc:
        status.update(stage='failed',error=str(exc));raise
    finally:
        if out is not None and out.is_dir():
            (out/'checkpoint_io_audit.json').write_text(json.dumps(status,indent=2),encoding='utf-8')
