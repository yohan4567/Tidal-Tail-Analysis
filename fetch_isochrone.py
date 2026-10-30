"""Request a documented, single-age Gaia EDR3/DR3 PARSEC isochrone from official CMD."""
import argparse
from html.parser import HTMLParser
import json
from pathlib import Path
import re
import urllib.parse
import urllib.request

from network import open_url
from provenance import sha256_file

ENDPOINT='https://stev.oapd.inaf.it/cgi-bin/cmd_3.8'


class Defaults(HTMLParser):
    def __init__(self):
        super().__init__(); self.values={}; self.select=None

    def handle_starttag(self,tag,attrs):
        a=dict(attrs)
        if tag=='select': self.select=a.get('name')
        if tag=='option' and self.select and 'selected' in a:
            self.values[self.select]=a['value']
        if tag=='input' and a.get('name') and a['name'] not in ('.cgifields','reset_form'):
            kind=a.get('type')
            if kind not in ('radio','checkbox') or 'checked' in a:
                self.values[a['name']]=a.get('value','')

    def handle_endtag(self,tag):
        if tag=='select': self.select=None


def validate_returned_model(table,age_years,mh):
    """Do not accept CMD boundary clipping as the requested composition or age."""
    import numpy as np
    for column,expected,tolerance in [('MH',mh,1e-4),('logAge',np.log10(age_years),1e-5)]:
        if column not in table or table[column].nunique(dropna=False)!=1:
            raise ValueError('Official isochrone lacks a single '+column)
        actual=float(table[column].iloc[0])
        if not np.isfinite(actual) or not np.isclose(actual,expected,atol=tolerance,rtol=0):
            raise ValueError(f'CMD returned {column}={actual}, requested {expected}; model-grid clipping/mismatch must be resolved, not silently relabelled')


def fetch(age_years,mh,out,*,av=0.,eta_reimers=.2,source):
    import numpy as np
    from population_model import read_isochrone
    if not np.isfinite([age_years,mh,av,eta_reimers]).all() or age_years<=0 or av<0 or not 0<=eta_reimers<=1 or not source.strip():
        raise ValueError('Finite age/[M/H]/Av, positive age, nonnegative Av and scientific source required')
    out=Path(out); meta_path=out.with_suffix('.json')
    if out.exists() or meta_path.exists(): raise ValueError('Use fresh output paths')
    with open_url(ENDPOINT) as response: html=response.read().decode()
    form=Defaults(); form.feed(html); params=form.values
    params.update(photsys_file='YBC_tab_mag_odfnew/tab_mag_gaiaEDR3.dat',
        imf_file='tab_imf/imf_kroupa_orig.dat',isoc_isagelog='0',
        isoc_agelow=str(age_years),isoc_ageupp=str(age_years),isoc_dage='0',
        isoc_ismetlog='1',isoc_metlow=str(mh),isoc_metupp=str(mh),isoc_dmet='0',
        extinction_av=str(av),eta_reimers=str(eta_reimers),output_kind='0',submit_form='Submit')
    body=urllib.parse.urlencode(params).encode()
    request=urllib.request.Request(ENDPOINT,data=body)
    with open_url(request,timeout=180) as response: result=response.read().decode()
    links=re.findall(r'href\s*=\s*[\"\']?([^\"\'\s>]*output\d+\.dat)',result,re.IGNORECASE)
    if not links: raise RuntimeError('CMD did not return an isochrone; response starts: '+re.sub('<[^>]+>',' ',result)[-1200:])
    url=urllib.parse.urljoin(ENDPOINT,links[0])
    if urllib.parse.urlparse(url).hostname!='stev.oapd.inaf.it': raise ValueError('Unexpected CMD output host')
    with open_url(url,timeout=120) as response: content=response.read()
    out.parent.mkdir(parents=True,exist_ok=True); out.write_bytes(content)
    table=read_isochrone(out)
    if len(table)<30: raise ValueError('Official isochrone unexpectedly short')
    validate_returned_model(table,age_years,mh)
    meta=dict(endpoint=ENDPOINT,download_url=url,request_parameters=params,sha256=sha256_file(out),
        age_years=age_years,metallicity_mh=mh,av=av,eta_reimers=eta_reimers,source_description=source,
        photometric_system='Gaia EDR3 passbands (same photometric system as DR3)',
        abundance_assumption='Scaled-solar [M/H]; not automatically spectroscopic [Fe/H] for alpha-enhanced stars',
        synthetic_calibration=False,n_nodes=len(table))
    meta_path.write_text(json.dumps(meta,indent=2),encoding='utf-8')
    return meta


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--age-years',type=float,required=True); p.add_argument('--mh',type=float,required=True)
    p.add_argument('--av',type=float,default=0.); p.add_argument('--source',required=True)
    p.add_argument('--eta-reimers',type=float,default=.2,help='RGB stellar-wind prescription; not the cluster escape rate')
    p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    print(json.dumps(fetch(a.age_years,a.mh,a.out,av=a.av,eta_reimers=a.eta_reimers,source=a.source),indent=2))
