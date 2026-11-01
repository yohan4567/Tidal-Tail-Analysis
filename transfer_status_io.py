"""Retry transient Windows status-file locks without changing scientific data."""
import json,os,time
from pathlib import Path

def write_json(path,value,*,status_file=False,attempts=25,delay=.2):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    text=json.dumps(value,indent=2)
    tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    tmp.write_text(text,encoding='utf-8')
    for attempt in range(attempts):
        try:tmp.replace(path);return
        except PermissionError:
            if attempt+1<attempts:time.sleep(delay)
    if not status_file:raise PermissionError('Could not replace locked JSON file: '+str(path))
    # Readers catch incomplete JSON and retry. Only mutable progress records
    # may use this fallback; scientific inputs/manifests remain strict.
    path.write_text(text,encoding='utf-8')
    tmp.unlink(missing_ok=True)
