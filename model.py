"""Backbone elastic-network benchmark; distances Angstrom, angles radians."""
from pathlib import Path
import numpy as np
from scipy.linalg import eigh
from scipy.spatial.distance import pdist, squareform

ROOT = Path(__file__).resolve().parent

def read_backbone(path, chain='A'):
    residues = {}
    for line in Path(path).read_text().splitlines():
        if not line.startswith('ATOM  ') or line[21] != chain or line[16] not in (' ', 'A'):
            continue
        atom = line[12:16].strip()
        if atom not in ('N', 'CA', 'C'):
            continue
        key = (int(line[22:26]), line[26])
        entry = residues.setdefault(key, {'name': line[17:20], 'atoms': {}})
        if atom in entry['atoms']:
            raise ValueError(f'Duplicate atom {key} {atom}')
        entry['atoms'][atom] = [float(line[i:i+8]) for i in (30, 38, 46)]
    keys = sorted(residues)
    if any(set(residues[k]['atoms']) != {'N', 'CA', 'C'} for k in keys):
        raise ValueError('Incomplete backbone')
    xyz = np.array([residues[k]['atoms'][a] for k in keys for a in ('N', 'CA', 'C')])
    return xyz, keys, [residues[k]['name'] for k in keys]

def load_pair():
    x, keys, names = read_backbone(ROOT/'data/1BA2.pdb')
    y, other, other_names = read_backbone(ROOT/'data/2DRI.pdb')
    if keys != other or any(b[0] != a[0]+1 for a,b in zip(keys,keys[1:])):
        raise ValueError('Residue correspondence/continuity requires manual preparation')
    bonds = np.linalg.norm(np.diff(x, axis=0), axis=1)
    if not np.all((bonds > 1.0) & (bonds < 1.8)):
        raise ValueError('Backbone break or invalid geometry')
    x -= x.mean(axis=0)
    y = align(y, x)[0]
    return x, y, names, [{'residue': k[0], 'start': a, 'target': b}
                          for k,a,b in zip(keys,names,other_names) if a != b]

def align(x, target):
    """Unweighted Kabsch alignment over CA atoms; row-vector convention."""
    ca = np.arange(1,len(x),3)
    xc, yc = x[ca].mean(0), target[ca].mean(0)
    u, _, vt = np.linalg.svd((x[ca]-xc).T @ (target[ca]-yc))
    correction = np.eye(3); correction[-1,-1] = np.linalg.det(u @ vt)
    rotation = u @ correction @ vt
    return (x-xc) @ rotation + yc, rotation

def objective(x, target):
    fitted, rotation = align(x, target)
    diff = fitted[1::3] - target[1::3]
    loss = np.mean(np.sum(diff*diff,axis=1))
    grad = np.zeros_like(x)
    grad[1::3] = (2/len(diff))*diff @ rotation.T
    return loss, grad

def rmsd(x,y):
    return np.sqrt(objective(x,y)[0])

def hessian(x, cutoff=8.0):
    n=len(x); h=np.zeros((3*n,3*n))
    ii,jj=np.where(np.triu(squareform(pdist(x)) < cutoff,1))
    for i,j in zip(ii,jj):
        delta=x[i]-x[j]; u=delta/np.linalg.norm(delta)
        # Stronger consecutive and angle-defining pairs; same in both methods.
        weight=10.0 if j-i==1 else (5.0 if j-i==2 else 1.0)
        block=weight*np.outer(u,u)
        a=slice(3*i,3*i+3); b=slice(3*j,3*j+3)
        h[a,a]+=block; h[b,b]+=block; h[a,b]-=block; h[b,a]-=block
    return h, len(ii)

def torsion_axes(names):
    axes=[]; labels=[]
    for r,name in enumerate(names):
        # First phi is an overall rigid rotation in this backbone representation.
        if r>0 and name!='PRO':
            axes.append(3*r+1); labels.append(f'phi_{r+1}')
        if r<len(names)-1:
            axes.append(3*r+2); labels.append(f'psi_{r+1}')
    return np.array(axes,dtype=int), labels

def reconstruct(x, axes, delta):
    out=x.copy()
    for i,angle in zip(axes,delta):
        if angle==0: continue
        p=out[i]; u=p-out[i-1]; u=u/np.linalg.norm(u)
        v=out[i+1:]-p
        c,s=np.cos(angle),np.sin(angle)
        out[i+1:]=p+c*v+s*np.cross(u,v)+(1-c)*np.outer(v@u,u)
    return out

def jacobian(x,axes):
    j=np.zeros((3*len(x),len(axes)))
    for k,i in enumerate(axes):
        u=x[i]-x[i-1];u/=np.linalg.norm(u)
        j[3*(i+1):,k]=np.cross(u,x[i+1:]-x[i]).ravel()
    return j

def project_rigid(x,j):
    xc=x-x.mean(0); basis=[]
    for u in np.eye(3):basis.append(np.tile(u,(len(x),1)).ravel())
    for u in np.eye(3):basis.append(np.cross(u,xc).ravel())
    q,_=np.linalg.qr(np.array(basis).T)
    return j-q@(q.T@j)

def torsion_gradient(x,axes,g):
    # Analytic suffix sums avoid building a full Jacobian on every iteration.
    suffix_g=np.cumsum(g[::-1],axis=0)[::-1]
    suffix_t=np.cumsum(np.cross(x,g)[::-1],axis=0)[::-1]
    u=x[axes]-x[axes-1];u/=np.linalg.norm(u,axis=1)[:,None]
    torque=suffix_t[axes+1]-np.cross(x[axes],suffix_g[axes+1])
    return np.sum(u*torque,axis=1)

def modes(x,names,method,k,cutoff=8.0):
    h,edges=hessian(x,cutoff)
    info={'atoms':len(x),'edges':edges,'cutoff_A':cutoff}
    axes,labels=torsion_axes(names)
    if method=='cartesian':
        vals,vec=eigh(h,subset_by_index=[0,min(len(h)-1,k+11)],driver='evr')
        pos=vals>max(1.0,np.max(vals))*1e-8
        info['zero_modes_in_low_spectrum']=int(np.sum(~pos))
        if info['zero_modes_in_low_spectrum']!=6:
            raise ValueError('Network does not have exactly six rigid-body zero modes')
        vals,basis=vals[pos][:k],vec[:,pos][:,:k]
        basis=basis.reshape(len(x),3,k)
        scale=np.sqrt(np.mean(np.sum(basis[1::3]**2,axis=1),axis=0))
        basis=basis/scale
    else:
        raw=jacobian(x,axes);j=project_rigid(x,raw)
        metric=j.T@j; hq=j.T@h@j
        mu,z=eigh(metric)
        keep=mu>mu[-1]*1e-10
        info.update(torsions=len(axes),retained_metric_rank=int(sum(keep)),metric_relative_cutoff=1e-10)
        w=z[:,keep]/np.sqrt(mu[keep])
        reduced=w.T@hq@w;reduced=(reduced+reduced.T)/2
        vals,v=eigh(reduced,subset_by_index=[0,k-1])
        if vals[0] <= 0:raise ValueError('Nonpositive internal mode')
        basis=w@v
        tangent=(j@basis).reshape(len(x),3,k)
        scale=np.sqrt(np.mean(np.sum(tangent[1::3]**2,axis=1),axis=0))
        basis=basis/scale
    info['eigenvalues']=vals.tolist()
    return basis,axes,info

def geometry(x,reference):
    lengths=np.linalg.norm(np.diff(x,axis=0),axis=1)
    ref_lengths=np.linalg.norm(np.diff(reference,axis=0),axis=1)
    def angles(z):
        a=z[:-2]-z[1:-1];b=z[2:]-z[1:-1]
        return np.arccos(np.clip(np.sum(a*b,axis=1)/(np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1)),-1,1))
    i,j=np.triu_indices(len(x),k=4)
    distance=np.linalg.norm(x[i]-x[j],axis=1)
    # Screening metric only; not an all-atom steric validation.
    bad=distance<2.0
    initial=np.linalg.norm(reference[i]-reference[j],axis=1)<2.0
    return {'max_backbone_bond_change_A':float(np.max(abs(lengths-ref_lengths))),
            'max_backbone_angle_change_deg':float(np.rad2deg(np.max(abs(angles(x)-angles(reference))))),
            'backbone_close_contacts_under_2A':int(bad.sum()),
            'new_backbone_close_contacts_under_2A':int(np.sum(bad&~initial))}
