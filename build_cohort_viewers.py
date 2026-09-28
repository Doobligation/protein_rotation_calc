"""Build offline per-protein motion viewers from saved benchmark amplitudes."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[name]='1'
import json,html,argparse
from pathlib import Path
import numpy as np
from model import ROOT,align
from cohort_model import prepare,basis,coordinates

def build(pairs,folder):
    folder=Path(folder);template=(ROOT/'viewer_template.html').read_text();cards=[]
    for index,pair in enumerate(pairs):
        names=['cartesian','torsion']
        stored=[folder/f'{index:02d}_forward_{name}_k10_seed0.json' for name in names]
        if not all(path.is_file() for path in stored):continue
        data=[json.loads(path.read_text()) for path in stored]
        if any(r.get('status')!='success' for r in data):continue
        try:
            x,y,names_res,segments,chains,_=prepare(pair)
            records=[]
            for r in data:
                B,axes,ends,_=basis(x,names_res,segments,chains,r['method'],r['k'])
                coef=np.array(r['amplitudes']);frames=[]
                for frac in np.linspace(0,1,31):
                    z=x+np.einsum('ijk,k->ij',B,frac*coef) if r['method']=='cartesian' else coordinates(x,axes,ends,chains,B@(frac*coef))
                    frames.append(align(z,y)[0][1::3].round(4).tolist())
                records.append({'method':r['method'],'k':r['k'],'seed':r['seed'],
                                'final_rmsd_A':r['final_rmsd_A'],'fractional_rmsd':r['fractional_rmsd'],
                                'total_seconds':r['total_seconds'],'setup_seconds':r['setup_seconds'],
                                'baseline_rmsd_A':r['baseline_rmsd_A'],'history':r['history'],'frames':frames,
                                'optimizer_success':r['status']=='success','optimizer_message':r.get('optimizer_message',''),
                                'geometry':{'max_backbone_bond_change_A':r['max_bond_change_A'],
                                            'max_backbone_angle_change_deg':r['max_angle_change_deg'],
                                            'new_backbone_close_contacts_under_2A':0}})
            payload={'start':x[1::3].round(4).tolist(),'target':y[1::3].round(4).tolist(),
                     'breaks':[begin for begin,_ in segments[1:]],'records':records}
            title=html.escape(pair['name']);text=template.replace('Ribose-binding protein',title)
            text=text.replace('1BA2 chain A → 2DRI chain A · 271 residues · N–Cα–C backbone model',
                              f"{html.escape(pair['a'])} chains {html.escape(pair['chains_a'])} → {html.escape(pair['b'])} chains {html.escape(pair['chains_b'])} · {len(names_res)} matched residues · N–Cα–C backbone")
            text=text.replace('1BA2 contains D67R; 2DRI is wild type and ribose-bound. Waters and ribose are excluded from this fitting model.',
                              'Matched structures may differ in sequence and ligands. Waters, ligands and side chains are excluded. See REPORT.md for missing backbone segments.')
            text=text.replace('● 2DRI reference',f"● {html.escape(pair['b'])} reference")
            text=text.replace(' · ${g.new_backbone_close_contacts_under_2A} new close contacts.',
                              ' · close contacts not assessed in this viewer.')
            text=text.replace('Close contacts use a simple backbone-only 2 Å screening threshold and are not all-atom clash validation.',
                              'Close-contact counts are unavailable in this cohort view. Check the protein geometry, residue coverage and chain breaks in REPORT.md and summary.csv.')
            text=text.replace('Timing is from this execution environment; three seeds sample fitting variability, not an extensive hardware performance study.',
                              'Timing is from one run on this machine. A single seed does not measure fitting variability.')
            text=text.replace('__DATA__',json.dumps(payload,separators=(',',':')))
            output=folder/f'protein_{index+1:02d}.html';output.write_text(text)
            cards.append(f'<li><a href="{output.name}">{index+1:02d}. {title}</a> — {len(names_res)} matched residues; {r["internal_chain_breaks"]} internal chain breaks</li>')
            print(output.name,flush=True)
        except Exception as error:
            print('VIEWER FAILED',pair['name'],type(error).__name__,str(error),flush=True)
    index_html='<!doctype html><html><meta charset="utf-8"><title>Protein benchmark viewers</title><style>body{font:16px system-ui;max-width:850px;margin:40px auto;background:#101821;color:#e7edf3}li{margin:12px}a{color:#81d8ce}</style><h1>Protein benchmark viewers</h1><p><a href="REPORT.md">Measured results and limitations</a></p><ol>'+''.join(cards)+'</ol></html>'
    (folder/'index.html').write_text(index_html)
    print(f'Viewers built: {len(cards)}/{len(pairs)}')
if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',default=str(ROOT/'pairs.json'))
    parser.add_argument('--output',default=str(ROOT/'results_13'))
    arguments=parser.parse_args()
    build(json.loads(Path(arguments.config).read_text()),Path(arguments.output))
