"""Numerical checks against independent finite differences and invariants."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np
from model import *
x,y,names,mismatch=load_pair();axes,_=torsion_axes(names)
assert len(x)==813 and len(names)==271 and len(mismatch)==1
assert mismatch[0]=={'residue':67,'start':'ARG','target':'ASP'}
assert rmsd(x,x)<1e-12
rng=np.random.default_rng(42)
rot,_=np.linalg.qr(rng.normal(size=(3,3)))
if np.linalg.det(rot)<0:rot[:,0]*=-1
assert rmsd(x@rot+np.array([3,-2,1]),x)<1e-10
assert np.allclose(reconstruct(x,axes,np.zeros(len(axes))),x)
q=rng.normal(0,.005,len(axes));z=reconstruct(x,axes,q)
g=geometry(z,x)
assert g['max_backbone_bond_change_A']<1e-10
assert g['max_backbone_angle_change_deg']<1e-8
loss,gx=objective(z,y);analytic=torsion_gradient(z,axes,gx)
errors=[]
for k in [0,10,len(axes)//2,len(axes)-1]:
 step=np.zeros(len(axes));step[k]=1e-6
 numerical=(objective(reconstruct(x,axes,q+step),y)[0]-objective(reconstruct(x,axes,q-step),y)[0])/2e-6
 errors.append(abs(numerical-analytic[k])/max(1.,abs(numerical)))
assert max(errors)<1e-5,errors
h,edges=hessian(x);assert np.allclose(h,h.T)
rigid=np.tile([1.,0.,0.],len(x));assert np.linalg.norm(h@rigid)<1e-9
v=rng.normal(size=x.shape);v/=np.linalg.norm(v)
# Independent finite difference of the nonlinear spring energy.
from scipy.spatial.distance import pdist,squareform
D=squareform(pdist(x));i,j=np.where(np.triu(D<8,1));w=np.where(j-i==1,10.,np.where(j-i==2,5.,1.))
def energy(z):return .5*np.sum(w*(np.linalg.norm(z[i]-z[j],axis=1)-D[i,j])**2)
step=1e-3
curvature=(energy(x+step*v)+energy(x-step*v))/step**2
assert abs(curvature-v.ravel()@h@v.ravel())<1e-5
print('PASS: correspondence, rigid alignment, zero-angle reconstruction, bond/angle preservation, torsion gradient, Hessian symmetry/nullspace/curvature.')
print('Initial CA RMSD:',rmsd(x,y),'A; torsions:',len(axes),'edges:',edges,'gradient relative errors:',errors)
