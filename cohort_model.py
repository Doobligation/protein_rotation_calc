"""Shared atom-matched backbone model for the paper's thirteen structure pairs."""
from collections import OrderedDict
from difflib import SequenceMatcher
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix, eye
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigh
from model import ROOT,align,rmsd,project_rigid

AA3 = dict(zip(
    'ALA ARG ASN ASP CYS GLU GLN GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL'.split(),
    'ARNDCEQGHILKMFPSTWYV'))

def structure_path(value):
    candidate = Path(value)
    if candidate.is_file():
        return candidate
    candidate = ROOT/'data'/f'{value}.pdb'
    if candidate.is_file():
        return candidate
    raise FileNotFoundError(f'PDB not found: {value}. Pass a local .pdb path or bundle a PDB ID in data/.')

def match_residues(a,b):
    """Prefer author residue keys; align sequences when numbering differs."""
    common=sorted(set(a)&set(b))
    identity = (sum(a[key]['name']==b[key]['name'] for key in common)/len(common)) if common else 0
    if len(common)>=.75*min(len(a),len(b)) and identity>=.8:
        return [(key,key) for key in common], 'author residue keys'
    ak=list(a);bk=list(b)
    sa=''.join(AA3.get(a[key]['name'],'X') for key in ak)
    sb=''.join(AA3.get(b[key]['name'],'X') for key in bk)
    pairs=[]
    for block in SequenceMatcher(None,sa,sb,autojunk=False).get_matching_blocks():
        pairs.extend((ak[block.a+i],bk[block.b+i]) for i in range(block.size))
    if len(pairs)<.75*min(len(a),len(b)):
        raise ValueError('Too few matching residues. These may be different proteins or require a curated sequence alignment.')
    return pairs, 'exact sequence blocks (different residue numbering)'

def residues(path,chain):
    out=OrderedDict()
    for line in Path(path).read_text().splitlines():
        if not line.startswith('ATOM  ') or line[21]!=chain or line[16] not in (' ','A'):
            continue
        atom=line[12:16].strip()
        if atom not in ('N','CA','C'):continue
        key=(int(line[22:26]),line[26])
        rec=out.setdefault(key,{'name':line[17:20],'atoms':{}})
        rec['atoms'][atom]=np.array([float(line[i:i+8]) for i in (30,38,46)])
    return {key:value for key,value in out.items() if set(value['atoms'])=={'N','CA','C'}}

def prepare(pair, return_metadata=False):
    ca=pair['chains_a'];cb=pair['chains_b']
    if ca=='auto':ca='ABCD'
    if cb=='auto':cb='ABCD'
    if len(ca)!=len(cb):raise ValueError('Chain count mismatch')
    X=[];Y=[];names=[];segments=[];chain_ranges=[];matches=[];match_details=[]
    for cha,chb in zip(ca,cb):
        a=residues(structure_path(pair['a']),cha);b=residues(structure_path(pair['b']),chb)
        if not a or not b:
            raise ValueError(f'Chain {cha}/{chb} has no complete N, CA, C residues')
        keys,match_type=match_residues(a,b)
        match_details.append({'chains':f'{cha}/{chb}','method':match_type,'matched':len(keys),
                              'input_a':len(a),'input_b':len(b)})
        start=len(names);last_x=None
        for key_a,key_b in keys:
            xa=np.array([a[key_a]['atoms'][t] for t in ('N','CA','C')]);yb=np.array([b[key_b]['atoms'][t] for t in ('N','CA','C')])
            if last_x is None or np.linalg.norm(xa[0]-last_x[-1])>1.8:
                segments.append([len(names),len(names)+1])
            else:segments[-1][1]+=1
            X.extend(xa);Y.extend(yb);names.append(a[key_a]['name']);matches.append((cha,chb,key_a[0],a[key_a]['name'],b[key_b]['name']))
            last_x=xa
        chain_ranges.append([start,len(names)])
    if len(names)<15:
        raise ValueError('At least 15 matched residues are required for this benchmark')
    if len(X)>3500:
        raise ValueError('This implementation limits matching backbones to 3500 atoms due to its dense torsion metric')
    x=np.array(X);y=np.array(Y);x-=x[1::3].mean(0);y=align(y,x)[0]
    if return_metadata:
        return x,y,names,segments,chain_ranges,matches,match_details
    return x,y,names,segments,chain_ranges,matches

def axes_for(names,segments):
    axes=[];ends=[]
    for beg,end in segments:
        for r in range(beg,end):
            if r>beg and names[r]!='PRO':axes.append(3*r+1);ends.append(3*end)
            if r<end-1:axes.append(3*r+2);ends.append(3*end)
    return np.array(axes,int),np.array(ends,int)

def coordinates(x,axes,ends,chain_ranges,q):
    out=x.copy();t=len(axes)
    for i,end,angle in zip(axes,ends,q[:t]):
        if not angle:continue
        u=out[i]-out[i-1];u/=np.linalg.norm(u)
        v=out[i+1:end]-out[i];c,s=np.cos(angle),np.sin(angle)
        out[i+1:end]=out[i]+c*v+s*np.cross(u,v)+(1-c)*np.outer(v@u,u)
    for c,(a,b) in enumerate(chain_ranges[1:]):
        start=3*a;end=3*b;angles=q[t+6*c:t+6*c+3];shift=q[t+6*c+3:t+6*c+6]
        if not np.any(angles) and not np.any(shift):continue
        z=out[start:end];center=z.mean(0);w=z-center
        for axis,theta in zip(np.eye(3),angles):
            if not theta:continue
            cosine,sine=np.cos(theta),np.sin(theta)
            w=cosine*w+sine*np.cross(axis,w)+(1-cosine)*np.outer(w@axis,axis)
        out[start:end]=center+w+shift
    return out

def jacobian(x,axes,ends,chain_ranges):
    n=len(x);d=len(axes)+6*(len(chain_ranges)-1);j=np.zeros((3*n,d))
    for k,(i,end) in enumerate(zip(axes,ends)):
        u=x[i]-x[i-1];u/=np.linalg.norm(u)
        j[3*(i+1):3*end,k]=np.cross(u,x[i+1:end]-x[i]).ravel()
    for c,(a,b) in enumerate(chain_ranges[1:]):
        first,last=3*a,3*b;center=x[first:last].mean(0);k=len(axes)+6*c
        for s,u in enumerate(np.eye(3)):
            j[3*first:3*last,k+s]=np.cross(u,x[first:last]-center).ravel()
            j[3*first:3*last,k+3+s]=np.tile(u,(last-first,1)).ravel()
    return project_rigid(x,j)

def gradient(x,axes,ends,chains,q,grad):
    values=np.zeros(len(q));sg=np.cumsum(grad[::-1],axis=0)[::-1];st=np.cumsum(np.cross(x,grad)[::-1],axis=0)[::-1]
    sg=np.vstack([sg,np.zeros((1,3))]);st=np.vstack([st,np.zeros((1,3))])
    for k,(i,end) in enumerate(zip(axes,ends)):
        u=x[i]-x[i-1];u/=np.linalg.norm(u)
        values[k]=u@(st[i+1]-st[end]-np.cross(x[i],sg[i+1]-sg[end]))
    for c,(a,b) in enumerate(chains[1:]):
        first,last=3*a,3*b;z=x[first:last];g=grad[first:last];center=z.mean(0)
        torque=np.cross(z-center,g).sum(0);k=len(axes)+6*c
        rx,ry,rz=q[k:k+3];cx,sx=np.cos(rx),np.sin(rx);cy,sy=np.cos(ry),np.sin(ry);cz,sz=np.cos(rz),np.sin(rz)
        Rz=np.array([[cz,-sz,0],[sz,cz,0],[0,0,1]])
        Ry=np.array([[cy,0,sy],[0,1,0],[-sy,0,cy]])
        # Effective final axes of sequential x, y, z rotations.
        values[k:k+3]=[torque@(Rz@Ry@np.array([1,0,0])),torque@(Rz@np.array([0,1,0])),torque[2]]
        values[k+3:k+6]=g.sum(axis=0)
    return values

def network(x,segments,cutoff=8.):
    n=len(x);p=np.array(list(cKDTree(x).query_pairs(cutoff)),int)
    if not len(p):raise ValueError('No contacts')
    i,j=p[:,0],p[:,1];seg=np.empty(n,int)
    for s,(a,b) in enumerate(segments):seg[3*a:3*b]=s
    sep=j-i;same=seg[i]==seg[j]
    w=np.where(same&(sep==1),10.,np.where(same&(sep==2),5.,1.))
    dv=x[i]-x[j];u=dv/np.linalg.norm(dv,axis=1)[:,None]
    blocks=w[:,None,None]*u[:,:,None]*u[:,None,:]
    a=np.arange(3);r=(3*i[:,None]+a).reshape(-1);c=(3*j[:,None]+a).reshape(-1)
    rr=[];cc=[];vv=[]
    for sgn,row,col in [(1,i,i),(1,j,j),(-1,i,j),(-1,j,i)]:
        rr.append((3*row[:,None,None]+a[None,:,None]+np.zeros((len(i),1,3),int)).ravel())
        cc.append((3*col[:,None,None]+a[None,None,:]+np.zeros((len(i),3,1),int)).ravel())
        vv.append((sgn*blocks).ravel())
    H=coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(3*n,3*n)).tocsr()
    return H,len(p)

def basis(x,names,segments,chains,method,k):
    h,edges=network(x,segments);axes,ends=axes_for(names,segments)
    if method=='cartesian':
        # Small-magnitude modes include six rigid-body zero modes; sparse eigensolve.
        # Shift/invert separates the exact rigid-body nullspace from very soft physical modes.
        # Adapt if disconnected fragments introduce additional zero modes.
        regularization=1e-6
        shifted=h+regularization*eye(h.shape[0],format='csr')
        count=k+8
        while True:
            vals,vec=eigsh(shifted,k=count,sigma=0.,which='LM',tol=1e-6,maxiter=3000)
            order=np.argsort(vals);vals=vals[order]-regularization;vec=vec[:,order]
            pos=vals>1e-8
            if pos.sum()>=k:break
            count=min(h.shape[0]-2,count+max(8,k))
            if count<=k or count>=h.shape[0]-2:raise ValueError('Insufficient positive Cartesian modes')
        vals=vals[pos][:k];B=vec[:,pos][:,:k].reshape(len(x),3,k)
        scale=np.sqrt(np.mean(np.sum(B[1::3]**2,axis=1),axis=0))
        return B/scale,axes,ends,{'edge_count':edges,'eigenvalues':vals.tolist(),'zero_modes':int(sum(~pos))}
    J=jacobian(x,axes,ends,chains);G=J.T@J
    mu,z=eigh(G,check_finite=False,overwrite_a=True);keep=mu>mu[-1]*1e-10
    W=z[:,keep]/np.sqrt(mu[keep]);JW=J@W
    Hq=JW.T@(h@JW);del J,JW
    vals,v=eigh(Hq,subset_by_index=(0,k-1),check_finite=False,overwrite_a=True)
    B=W@v
    # Need unprojected tangent for scaling (global rotations removed before fitting).
    tangent=jacobian(x,axes,ends,chains)@B
    tangent=tangent.reshape(len(x),3,k);scale=np.sqrt(np.mean(np.sum(tangent[1::3]**2,axis=1),axis=0))
    return B/scale,axes,ends,{'edge_count':edges,'torsions':len(axes),'rigid_chain_dofs':6*(len(chains)-1),'rank':int(keep.sum()),'eigenvalues':vals.tolist()}
