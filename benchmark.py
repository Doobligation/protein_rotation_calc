"""Run isolated, single-threaded Cartesian/torsion backbone benchmarks."""
import os
# Set before importing NumPy; subprocesses use the same declared thread count.
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key]='1'
import argparse,json,time,sys,subprocess,platform,hashlib,csv
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import scipy
from model import ROOT,load_pair,modes,reconstruct,objective,torsion_gradient,rmsd,geometry,align

def peak_memory_mb():
    try:
        import resource
        rss=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        return rss/(1024**2 if sys.platform=='darwin' else 1024)
    except ImportError:return None

def worker(args):
    start=time.perf_counter()
    x,y,names,mismatches=load_pair()
    basis,axes,info=modes(x,names,args.method,args.k,args.cutoff)
    setup=time.perf_counter()-start
    baseline=rmsd(x,y);history=[];best=[float('inf'),None]
    rng=np.random.default_rng(args.seed)
    initial=np.zeros(args.k) if args.seed==0 else rng.normal(0,0.25,args.k)
    def evaluate(a):
        if args.method=='cartesian':z=x+np.einsum('ijk,k->ij',basis,a)
        else:z=reconstruct(x,axes,basis@a)
        loss,g=objective(z,y)
        if args.method=='cartesian':grad=np.einsum('ijk,ij->k',basis,g)
        else:grad=basis.T@torsion_gradient(z,axes,g)
        if loss<best[0]:best[:]=[float(loss),a.copy()]
        history.append([time.perf_counter()-start,float(np.sqrt(best[0]))])
        return loss,grad
    initial_rmsd=float(np.sqrt(evaluate(initial)[0]))
    fit_start=time.perf_counter()
    result=minimize(evaluate,initial,jac=True,method='L-BFGS-B',bounds=[(-8.,8.)]*args.k,
                    options={'maxiter':args.maxiter,'ftol':1e-10,'gtol':1e-6,'maxls':30})
    fit_seconds=time.perf_counter()-fit_start
    a=best[1]
    z=x+np.einsum('ijk,k->ij',basis,a) if args.method=='cartesian' else reconstruct(x,axes,basis@a)
    screening=geometry(z,x);final=rmsd(z,y)
    total=time.perf_counter()-start
    # Path is a mode-amplitude interpolation, never a physical time trajectory.
    display_start=time.perf_counter()
    frames=[]
    for fraction in np.linspace(0,1,31):
        frame=x+np.einsum('ijk,k->ij',basis,fraction*a) if args.method=='cartesian' else reconstruct(x,axes,basis@(fraction*a))
        frames.append(align(frame,y)[0][1::3].round(4).tolist())
    display_seconds=time.perf_counter()-display_start
    record={'method':args.method,'k':args.k,'seed':args.seed,'baseline_rmsd_A':baseline,
            'initial_search_rmsd_A':initial_rmsd,'final_rmsd_A':final,'fractional_rmsd':final/baseline,
            'setup_seconds':setup,'fit_seconds':fit_seconds,'total_seconds':total,
            'display_frame_generation_seconds_excluded':display_seconds,'peak_process_rss_MB':peak_memory_mb(),
            'optimizer_success':bool(result.success),'optimizer_message':str(result.message),
            'iterations':int(result.nit),'function_evaluations':len(history),
            'amplitudes':a.tolist(),'amplitude_bounds':[-8,8],
            'geometry':screening,'model':info,'sequence_mismatches':mismatches,'history':history,'frames':frames,
            'coordinates':z.tolist()}
    # ru_maxrss includes viewer-frame preparation; include explicit scope for readers.
    record['peak_memory_scope']='whole worker including frame generation; interpreter baseline included'
    output=ROOT/'results'/f'{args.method}_k{args.k}_seed{args.seed}.json'
    output.write_text(json.dumps(record))
    print(f'{args.method:9s} K={args.k:2d} seed={args.seed} RMSD={final:.4f} A fRMS={final/baseline:.4f} time={total:.2f}s success={result.success}',flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--worker',action='store_true');parser.add_argument('--method',choices=['cartesian','torsion'])
    parser.add_argument('--k',type=int,default=10);parser.add_argument('--seed',type=int,default=0)
    parser.add_argument('--modes',type=int,nargs='+',default=[1,2,5,10]);parser.add_argument('--seeds',type=int,nargs='+',default=[0,1,2])
    parser.add_argument('--maxiter',type=int,default=150);parser.add_argument('--cutoff',type=float,default=8.)
    args=parser.parse_args();(ROOT/'results').mkdir(exist_ok=True)
    if args.worker:return worker(args)
    outputs=[]
    for method in ['cartesian','torsion']:
        for k in args.modes:
            for seed in args.seeds:
                subprocess.run([sys.executable,str(Path(__file__).resolve()),'--worker','--method',method,'--k',str(k),
                                '--seed',str(seed),'--maxiter',str(args.maxiter),'--cutoff',str(args.cutoff)],check=True)
                outputs.append(json.loads((ROOT/'results'/f'{method}_k{k}_seed{seed}.json').read_text()))
    manifest={'platform':platform.platform(),'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
              'threads':1,'parameters':vars(args),'inputs':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/'data').glob('*.pdb')},
              'runs':[{'method':r['method'],'k':r['k'],'seed':r['seed']} for r in outputs]}
    (ROOT/'results/manifest.json').write_text(json.dumps(manifest,indent=2))
    fields=['method','k','seed','final_rmsd_A','fractional_rmsd','setup_seconds','fit_seconds','total_seconds','peak_process_rss_MB','optimizer_success','iterations','function_evaluations']
    with (ROOT/'results/summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        writer.writerows({key:r[key] for key in fields} for r in outputs)
    from viewer import build
    build(outputs)

if __name__=='__main__':main()
