"""Measured backbone benchmark on all 13 named protein structure pairs."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS'):os.environ[key]='1'
import sys,json,time,csv,hashlib,platform,subprocess,argparse
from pathlib import Path
import numpy as np
from scipy.optimize import minimize
import scipy
from model import ROOT,objective,rmsd
from cohort_model import prepare,basis,coordinates,gradient

PAIRS=json.loads((ROOT/'pairs.json').read_text())
OUT=ROOT/'results_13'

def configure(args):
    global PAIRS, OUT
    if args.config and args.start:
        raise ValueError('Use --config or --start/--target, not both')
    if args.config:
        PAIRS=json.loads(Path(args.config).read_text())
        OUT=Path(args.output or ROOT/'results_custom')
    elif args.start:
        if not args.target or not args.start_chains or not args.target_chains:
            raise ValueError('--start requires --target, --start-chains and --target-chains')
        PAIRS=[{'name':args.name or 'Custom protein pair','a':args.start,'b':args.target,
                'chains_a':args.start_chains,'chains_b':args.target_chains}]
        OUT=Path(args.output or ROOT/'results_custom')
    elif args.output:
        OUT=Path(args.output)
    for pair in PAIRS:
        if not all(k in pair for k in ('name','a','b','chains_a','chains_b')):
            raise ValueError('Each pair requires name, a, b, chains_a and chains_b')

def memory():
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024**2 if sys.platform=='darwin' else 1024)
    except ImportError:return None

def work(args):
    pair=PAIRS[args.pair]
    if args.reverse:
        pair={**pair,'a':pair['b'],'b':pair['a'],'chains_a':pair['chains_b'],'chains_b':pair['chains_a']}
    method=args.method;t0=time.perf_counter()
    record={'pair':PAIRS[args.pair]['name'],'a':pair['a'],'b':pair['b'],'chains_a':pair['chains_a'],'chains_b':pair['chains_b'],
            'method':method,'k':args.k,'seed':args.seed,'reverse':args.reverse,'status':'failed'}
    output=OUT/f"{args.pair:02d}_{'reverse' if args.reverse else 'forward'}_{method}_k{args.k}_seed{args.seed}.json"
    try:
        x,y,names,segments,chains,matches,match_details=prepare(pair,return_metadata=True)
        record.update(residues=len(names),atoms=len(x),segments=len(segments),chains=len(chains),
                      matching=match_details,internal_chain_breaks=len(segments)-len(chains),
                      sequence_mismatches=[list(m) for m in matches if m[3]!=m[4]],
                      baseline_rmsd_A=rmsd(x,y))
        B,axes,ends,details=basis(x,names,segments,chains,method,args.k)
        setup=time.perf_counter()-t0
        Bcart=B if method=='cartesian' else None
        rng=np.random.default_rng(args.seed)
        init=np.zeros(args.k) if args.seed==0 else rng.normal(0,.25,args.k)
        best=[float('inf'),None];evaluations=0;history=[]
        def calc(a):
            nonlocal evaluations
            z=x+np.einsum('ijk,k->ij',Bcart,a) if Bcart is not None else coordinates(x,axes,ends,chains,B@a)
            loss,g=objective(z,y)
            grad=np.einsum('ijk,ij->k',Bcart,g) if Bcart is not None else B.T@gradient(z,axes,ends,chains,B@a,g)
            if loss<best[0]:best[:]=[float(loss),a.copy()]
            evaluations+=1;history.append([time.perf_counter()-t0,float(np.sqrt(best[0]))])
            return loss,grad
        opt=minimize(calc,init,jac=True,method='L-BFGS-B',bounds=[(-8,8)]*args.k,
                     options={'maxiter':args.maxiter,'ftol':1e-10,'gtol':1e-6,'maxls':30})
        a=best[1];z=x+np.einsum('ijk,k->ij',Bcart,a) if Bcart is not None else coordinates(x,axes,ends,chains,B@a)
        bond_d=[];angle_d=[]
        for beg,end in segments:
            sl=slice(3*beg,3*end)
            v0=x[sl];v=z[sl]
            bond_d.extend(abs(np.linalg.norm(np.diff(v,axis=0),axis=1)-np.linalg.norm(np.diff(v0,axis=0),axis=1)))
            if len(v)<3:continue
            def angles(w):
                u=w[:-2]-w[1:-1];v=w[2:]-w[1:-1]
                return np.arccos(np.clip((u*v).sum(axis=1)/np.linalg.norm(u,axis=1)/np.linalg.norm(v,axis=1),-1,1))
            angle_d.extend(abs(angles(v)-angles(v0)))
        record.update(status='success' if opt.success else 'optimizer_stopped',
                      optimizer_message=str(opt.message),mode_metadata=details,setup_seconds=setup,
                      total_seconds=time.perf_counter()-t0,peak_rss_MB=memory(),
                      final_rmsd_A=rmsd(z,y),fractional_rmsd=rmsd(z,y)/record['baseline_rmsd_A'],
                      function_evaluations=evaluations,iterations=int(opt.nit),
                      max_bond_change_A=max(bond_d,default=0.),max_angle_change_deg=float(np.rad2deg(max(angle_d,default=0.))),
                      amplitudes=a.tolist(),history=history)
    except Exception as e:
        record.update(error_type=type(e).__name__,error=str(e),total_seconds=time.perf_counter()-t0)
    OUT.mkdir(exist_ok=True);output.write_text(json.dumps(record))
    print(f"{args.pair+1:02d} {'B→A' if args.reverse else 'A→B'} {method:9s} {record['status']:17s} {record.get('final_rmsd_A','-')} time={record['total_seconds']:.1f}s {record.get('error','')}",flush=True)
    if record['status']=='failed':sys.exit(2)

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--worker',action='store_true');ap.add_argument('--pair',type=int);ap.add_argument('--reverse',action='store_true')
    ap.add_argument('--method',choices=['cartesian','torsion']);ap.add_argument('--k',type=int,default=10)
    ap.add_argument('--maxiter',type=int,default=150);ap.add_argument('--seed',type=int,default=0)
    ap.add_argument('--both-directions',action='store_true');ap.add_argument('--pairs',type=int,nargs='+',default=list(range(13)))
    ap.add_argument('--config',help='JSON list of structure pairs')
    ap.add_argument('--start',help='PDB file path or bundled PDB ID')
    ap.add_argument('--target',help='PDB file path or bundled PDB ID')
    ap.add_argument('--start-chains',help='Matching chain IDs, e.g. A or AB')
    ap.add_argument('--target-chains',help='Corresponding target chain IDs')
    ap.add_argument('--name',help='Name for a custom pair')
    ap.add_argument('--output',help='Output directory')
    args=ap.parse_args()
    configure(args)
    if args.worker:return work(args)
    OUT.mkdir(exist_ok=True)
    if args.start or args.config:
        args.pairs=list(range(len(PAIRS)))
    runs=[]
    for i in args.pairs:
        for reverse in ([False,True] if args.both_directions else [False]):
            for method in ['cartesian','torsion']:
                cmd=[sys.executable,str(Path(__file__).resolve()),'--worker','--pair',str(i),'--method',method,'--k',str(args.k),'--maxiter',str(args.maxiter),'--seed',str(args.seed)]
                if reverse:cmd.append('--reverse')
                if args.config:cmd+=['--config',args.config]
                if args.start:cmd+=['--start',args.start,'--target',args.target,'--start-chains',args.start_chains,'--target-chains',args.target_chains]
                if args.name:cmd+=['--name',args.name]
                if args.output:cmd+=['--output',args.output]
                sub=subprocess.run(cmd)
                path=OUT/f"{i:02d}_{'reverse' if reverse else 'forward'}_{method}_k{args.k}_seed{args.seed}.json"
                runs.append(json.loads(path.read_text()) if path.exists() else {'pair':PAIRS[i]['name'],'method':method,'status':'failed','error':f'worker exited {sub.returncode}'})
                write_results(runs,args)
    if args.k==10 and args.seed==0:
        from cohort_report import report
        from build_cohort_viewers import build
        report(OUT,PAIRS)
        build(PAIRS,OUT)

def write_results(runs,args):
    fields=['pair','a','b','chains_a','chains_b','method','k','seed','reverse','status','residues','segments','baseline_rmsd_A','final_rmsd_A','fractional_rmsd','setup_seconds','total_seconds','peak_rss_MB','max_bond_change_A','max_angle_change_deg','function_evaluations','error']
    with (OUT/'summary.csv').open('w',newline='') as f:
        wr=csv.DictWriter(f,fields);wr.writeheader();wr.writerows({k:r.get(k,'') for k in fields} for r in runs)
    used={value for p in PAIRS for value in (p['a'],p['b'])}
    inputs={str(value):hashlib.sha256((Path(value) if Path(value).is_file() else ROOT/'data'/f'{value}.pdb').read_bytes()).hexdigest()
            for value in used}
    data={'parameters':vars(args),'platform':platform.platform(),'python':sys.version,'numpy':np.__version__,'scipy':scipy.__version__,
          'inputs':inputs,
          'statuses':{p['name']:[r['status'] for r in runs if r['pair']==p['name']] for p in PAIRS}}
    (OUT/'manifest.json').write_text(json.dumps(data,indent=2))
if __name__=='__main__':main()
