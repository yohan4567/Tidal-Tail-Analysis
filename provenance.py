"""Small streaming helpers for recording large input catalogues."""
import hashlib


def sha256_file(path):
    digest=hashlib.sha256()
    with open(path,'rb') as f:
        for chunk in iter(lambda:f.read(1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()
