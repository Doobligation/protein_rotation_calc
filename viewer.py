"""Generate an offline, dependency-free HTML viewer from measured results."""
import json
from model import ROOT,load_pair

def build(records):
    x,y,_,_=load_pair()
    payload={'start':x[1::3].round(4).tolist(),'target':y[1::3].round(4).tolist(),
             'records':[{k:v for k,v in r.items() if k!='coordinates'} for r in records]}
    template=(ROOT/'viewer_template.html').read_text()
    (ROOT/'results/viewer.html').write_text(template.replace('__DATA__',json.dumps(payload,separators=(',',':'))))
    print('Viewer: results/viewer.html',flush=True)

if __name__=='__main__':
    manifest=json.loads((ROOT/'results/manifest.json').read_text())
    records=[json.loads((ROOT/'results'/f"{r['method']}_k{r['k']}_seed{r['seed']}.json").read_text()) for r in manifest['runs']]
    build(records)
