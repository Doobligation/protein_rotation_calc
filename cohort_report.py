"""Summarize the actual cohort outputs, recording exclusions and limitations."""
from pathlib import Path
import csv,json,math
import numpy as np
from model import ROOT
from benchmark_13 import PAIRS

def report(folder=ROOT/'results_13',pairs=PAIRS):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    folder=Path(folder)
    records={}
    for f in folder.glob('??_*_k*_seed*.json'):
        r=json.loads(f.read_text())
        records[(r['pair'],r['reverse'],r['method'])]=r
    directions=[False,True] if any(k[1] for k in records) else [False]
    rows=[]
    for p in pairs:
        for reverse in directions:
            c=records.get((p['name'],reverse,'cartesian'))
            t=records.get((p['name'],reverse,'torsion'))
            rows.append((p,reverse,c,t))
    def fmt(r,k,precision=3):
        return f"{r[k]:.{precision}f}" if r and r.get('status')=='success' and r.get(k) is not None else '—'
    successful=[r for r in records.values() if r.get('status')=='success']
    failed=[r for r in records.values() if r.get('status')!='success']
    complete=[]
    for p in pairs:
        f=records.get((p['name'],False,'torsion'))
        if f and f.get('status')=='success' and not f['internal_chain_breaks']:
            complete.append(p['name'])
    paper_count=len(pairs)
    lines=['# Thirteen protein structure pairs: measured benchmark','',
           f"Attempted {len(records)} method × direction runs across {paper_count} structure pairs. "
           f"{len(successful)} converged, {len(failed)} did not converge or failed. "
           'These are N–Cα–C backbone fits with the target provided during optimization.','',
           '## Individual results','',
           '| Protein | Direction | Matched residues | Internal chain breaks | Initial RMSD (Å) | Cartesian fRMS | Torsion fRMS | Cartesian time (s) | Torsion time (s) |',
           '|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for p,rev,c,t in rows:
        r=c or t
        lines.append(f"| {p['name']} | {'B→A' if rev else 'A→B'} | {r.get('residues','—') if r else '—'} | {r.get('internal_chain_breaks','—') if r else '—'} | {fmt(r,'baseline_rmsd_A',2)} | {fmt(c,'fractional_rmsd')} | {fmt(t,'fractional_rmsd')} | {fmt(c,'total_seconds',2)} | {fmt(t,'total_seconds',2)} |")
    direct=[(p,rev,c,t) for p,rev,c,t in rows if c and t and c.get('status')=='success' and t.get('status')=='success']
    better=sum(t['fractional_rmsd']<c['fractional_rmsd'] for _,_,c,t in direct)
    lines+=['','![Measured comparison](comparison_13.png)','',
            f'Across {len(direct)} directly paired direction tests, the torsion implementation had lower fitted Cα fractional RMSD in {better}. '
            'The number describes these inputs and this prototype; it is not a statistical claim about all proteins.','',
            '## Scientific limits','',
            '- Pair matching uses author residue keys when they cover at least 75% of the smaller structure with at least 80% identical residues; otherwise it uses exact sequence blocks. Missing and unmatched residues are excluded from both structures. See each run JSON for the matching method and sequence differences.',
            '- A chain break is recorded when adjacent matched backbone segments have a C–N distance over 1.8 Å. No torsion is assigned across that break. These pairs cannot be compared to the paper as full-chain models. Multi-chain models include six relative rigid-body degrees of freedom per additional chain.',
            '- The paper used different atom and mode models and a different optimizer. Do not compare these numeric fRMS values to published all-atom fRMS as a reproduction.',
            '- Cartesian outputs can stretch bonds and change angles. The summary CSV records each model’s maximum backbone changes. A low RMSD is a structural fitting score and is not a physical conformer or a free-energy estimate.',
            '- Timings use one CPU thread and include PDB preparation, network, modes, fit and final backbone geometry check. Cache, process startup and optional viewer generation are excluded; these timings are environment-specific. Peak RSS includes Python and imported libraries.',
            '- One initialization seed and ten modes are used for each direction. Larger samples, seed repeats and cutoff sensitivity checks are needed before general performance claims.',
            '',f'Forward pairs without internal matched-chain breaks: {len(complete)}/{paper_count} ({", ".join(complete)}).',
            '', '## Reproduce','',
            'From the project directory with dependencies installed:','',
            '```bash','python benchmark_13.py --both-directions --k 10 --seed 0','python cohort_report.py','```','',
            'To run another protein pair, see the README and pass local PDB paths and corresponding chain IDs to `benchmark_13.py`.','']
    if failed:
        lines+=['','## Nonconverged or failed runs','']
        for r in failed:lines.append(f"- {r['pair']} {'B→A' if r['reverse'] else 'A→B'} {r['method']}: {r.get('error',r.get('optimizer_message','unknown'))}")
    (folder/'REPORT.md').write_text('\n'.join(lines))
    colors={'cartesian':'#277caa','torsion':'#8d51b2'}
    fig,(ax,ax2)=plt.subplots(1,2,figsize=(15,7),layout='constrained')
    data=[(p,records.get((p['name'],False,'cartesian')),records.get((p['name'],False,'torsion'))) for p in pairs]
    y=np.arange(len(data))
    for method,axis,shift in [('cartesian',ax,-.16),('torsion',ax,.16)]:
        vals=[r.get('fractional_rmsd',np.nan) if r and r.get('status')=='success' else np.nan for p,c,t in data for r in [c if method=='cartesian' else t]]
        axis.scatter(vals,y+shift,label=method.title(),color=colors[method],s=32)
    ax.set_yticks(y,[p['name'] for p,_,_ in data]);ax.set_xlabel('Fractional Cα RMSD (lower is better)');ax.axvline(1,ls=':',color='gray');ax.invert_yaxis();ax.grid(axis='x',alpha=.2);ax.legend()
    for method,shift in [('cartesian',-.16),('torsion',.16)]:
        vals=[r.get('total_seconds',np.nan) if r and r.get('status')=='success' else np.nan for p,c,t in data for r in [c if method=='cartesian' else t]]
        ax2.scatter(vals,y+shift,label=method.title(),color=colors[method],s=32)
    ax2.set_yticks(y,['']*len(data));ax2.set_xlabel('Measured computation time (s)');ax2.set_xscale('log');ax2.invert_yaxis();ax2.grid(axis='x',alpha=.2)
    fig.suptitle('13 protein pairs · A→B · ten modes · backbone fits',fontsize=15)
    fig.savefig(folder/'comparison_13.png',dpi=160);plt.close(fig)
    print(f'{len(successful)}/{len(records)} successful; torsion lower fRMS on {better}/{len(direct)} paired directions. Report: {folder/"REPORT.md"}')
if __name__=='__main__':report()
