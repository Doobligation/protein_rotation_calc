"""Create summary plots and a Markdown report from completed measured runs."""
import json,csv
import numpy as np
from model import ROOT

def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    manifest=json.loads((ROOT/'results/manifest.json').read_text())
    records=[json.loads((ROOT/'results'/f"{r['method']}_k{r['k']}_seed{r['seed']}.json").read_text()) for r in manifest['runs']]
    ks=sorted({r['k'] for r in records});rows=[]
    for method in ['cartesian','torsion']:
        for k in ks:
            group=[r for r in records if r['method']==method and r['k']==k]
            rows.append({'method':method,'modes':k,'median_rmsd_A':float(np.median([r['final_rmsd_A'] for r in group])),
                         'median_time_seconds':float(np.median([r['total_seconds'] for r in group])),
                         'min_time_seconds':min(r['total_seconds'] for r in group),'max_time_seconds':max(r['total_seconds'] for r in group),
                         'median_peak_rss_MB':float(np.median([r['peak_process_rss_MB'] for r in group if r['peak_process_rss_MB'] is not None])) if any(r['peak_process_rss_MB'] is not None for r in group) else None,
                         'max_bond_change_A':max(r['geometry']['max_backbone_bond_change_A'] for r in group),
                         'max_angle_change_deg':max(r['geometry']['max_backbone_angle_change_deg'] for r in group),
                         'max_new_close_contacts':max(r['geometry']['new_backbone_close_contacts_under_2A'] for r in group)})
    with (ROOT/'results/aggregate.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    thresholds=[]
    for r in records:
        for threshold in [2.,1.5,1.]:
            found=[t for t,e in r['history'] if e<=threshold]
            thresholds.append({'method':r['method'],'modes':r['k'],'seed':r['seed'],'rmsd_threshold_A':threshold,
                               'time_seconds':min(found) if found else '', 'status':'reached' if found else 'not reached',
                               'scope':'RMSD only; geometry not screened at every history point'})
    with (ROOT/'results/threshold_times.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(thresholds[0]));writer.writeheader();writer.writerows(thresholds)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained')
    for method,color in [('cartesian','#277caa'),('torsion','#8955b0')]:
        group=[r for r in rows if r['method']==method]
        axs[0].plot(ks,[r['median_rmsd_A'] for r in group],'o-',label=method.title(),color=color)
        vals=np.array([r['median_time_seconds'] for r in group]);low=np.array([r['min_time_seconds'] for r in group]);high=np.array([r['max_time_seconds'] for r in group])
        axs[1].errorbar(ks,vals,yerr=[vals-low,high-vals],fmt='o-',capsize=3,label=method.title(),color=color)
    axs[0].axhline(records[0]['baseline_rmsd_A'],ls=':',color='gray',label='Starting structure')
    axs[0].set(ylabel='Aligned Cα RMSD (Å)',xlabel='Number of modes',title='Fit accuracy (median across seeds)')
    axs[1].set(ylabel='Wall time (seconds)',xlabel='Number of modes',title='Compute time (median and range)')
    for ax in axs:ax.set_xticks(ks);ax.grid(alpha=.2);ax.legend(fontsize=8)
    fig.suptitle('Ribose-binding protein: backbone fitting benchmark',fontsize=14)
    fig.savefig(ROOT/'results/comparison.png',dpi=170);plt.close(fig)
    lines=['# Measured results: ribose-binding protein','',
           'Backbone-only N–Cα–C elastic-network fitting of 1BA2 chain A to 2DRI chain A. This is a target-guided structural benchmark, not molecular dynamics or a free-energy calculation.','',
           f"Starting aligned Cα RMSD: **{records[0]['baseline_rmsd_A']:.6f} Å**. Both chains contain 271 matched residues and 813 modeled atoms. The only residue identity mismatch is ARG67 versus ASP67.",'',
           '| Method | Modes | Median RMSD (Å) | Median compute time (s) | Time range (s) | Max bond change (Å) | Max angle change (degrees) |',
           '|---|---:|---:|---:|---:|---:|---:|']
    for r in rows:lines.append(f"| {r['method']} | {r['modes']} | {r['median_rmsd_A']:.4f} | {r['median_time_seconds']:.3f} | {r['min_time_seconds']:.3f}–{r['max_time_seconds']:.3f} | {r['max_bond_change_A']:.4g} | {r['max_angle_change_deg']:.4g} |")
    lines+=['','![Measured comparison](comparison.png)','','## Interpretation','',
            'The torsion model reaches lower Cα fitting error at each tested mode count in this run and preserves the starting backbone bond lengths and angles. The Cartesian fits distort local backbone geometry substantially. A low RMSD alone therefore does not establish a physically usable conformer.',
            '', 'These are measurements from one machine with three initialization seeds per configuration. The methods use the same target, atom selection, contact potential, L-BFGS-B implementation, iteration cap, and amplitude bounds, but they explore different constrained spaces. This measures the implemented pipelines; it does not isolate coordinate representation alone or prove that one algorithm is universally faster.',
            '', 'No all-atom relaxation is included. Neither side chains, hydrogens, solvent nor ribose are modeled. Geometry screens only concern the represented backbone; do not interpret zero backbone close contacts as an all-atom validation.',
            '', 'Peak RSS is the whole fresh worker-process maximum, including imports, mode construction, fitting and display-frame generation. It is not incremental algorithm memory, and is unavailable on Windows with the standard-library implementation. Compute time excludes display-frame generation and output writing. Parent process launch and disk cache effects are not included in compute time.',
            '', 'The chronological histories record the best objective value among optimizer evaluations, including trial points. threshold_times.csv gives RMSD-only threshold times; geometry was evaluated at the final result rather than every history point. Do not interpret these as times to a fully validated structure.',
            '', f"Optimizer success: {sum(r['optimizer_success'] for r in records)}/{len(records)} runs. Eigenmode calculations found six Cartesian rigid-body zero modes; torsion rank and spectrum are recorded per run.",
            '', '## Reproducibility','', 'See manifest.json for platform, dependency versions, single-thread settings, input SHA-256 hashes, seeds and run arguments. Run python check_model.py, python benchmark.py and python report.py from the project directory. Timings will vary on other computers.','',
            '## Files','', '- viewer.html: standalone interactive 3D Cα traces and measured convergence curves.', '- summary.csv / aggregate.csv: individual and grouped measurements.', '- threshold_times.csv: threshold success/failure with no fabricated times.', '- *_k*_seed*.json: full coordinates, display frames, geometry, modes metadata and histories.','']
    (ROOT/'results/RESULTS.md').write_text('\n'.join(lines))
    print('\n'.join(lines[:17]))
if __name__=='__main__':main()
