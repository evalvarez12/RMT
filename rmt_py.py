#RMT python

from pylab import *

#NOTA: En kron(a,b) a queda al final del estado 

SIGMA_X=array([[0,1],[1,0]])
SIGMA_Y=array([[0,-1j],[1j,0]])
SIGMA_Z=array([[1,0],[0,-1]])

def purity_C(state) :
  l=size(state)
  a=state[:l/2]
  b=state[l/2:]
  rho=array([[dot(a,a.conjugate()),dot(a,b.conjugate())],[dot(b,a.conjugate()),dot(b,b.conjugate())]])
  return trace(dot(rho,rho))

def purity_C2(state) :
  l=size(state)
  a=state[:l/2]
  a1=a[:l/4]
  a2=a[l/4:]
  b=state[l/2:]
  b1=b[:l/4]
  b2=b[l/4:]
  rho=array([[dot(a1,a1.conjugate()),dot(a1,a2.conjugate()),dot(a1,b1.conjugate()),dot(a1,b2.conjugate())],[dot(a2,a1.conjugate()),dot(a2,a2.conjugate()),dot(a2,b1.conjugate()),dot(a2,b2.conjugate())],[dot(b1,a1.conjugate()),dot(b1,a2.conjugate()),dot(b1,b1.conjugate()),dot(b1,b2.conjugate())],[dot(b2,a1.conjugate()),dot(b2,a2.conjugate()),dot(b2,b1.conjugate()),dot(b2,b2.conjugate())]])
  return trace(dot(rho,rho))

def concurrence(state) :
  l=size(state)
  a=state[:l/2]
  a1=a[:l/4]
  a2=a[l/4:]
  b=state[l/2:]
  b1=b[:l/4]
  b2=b[l/4:]
  rho=array([[dot(a1,a1.conjugate()),dot(a1,a2.conjugate()),dot(a1,b1.conjugate()),dot(a1,b2.conjugate())],[dot(a2,a1.conjugate()),dot(a2,a2.conjugate()),dot(a2,b1.conjugate()),dot(a2,b2.conjugate())],[dot(b1,a1.conjugate()),dot(b1,a2.conjugate()),dot(b1,b1.conjugate()),dot(b1,b2.conjugate())],[dot(b2,a1.conjugate()),dot(b2,a2.conjugate()),dot(b2,b1.conjugate()),dot(b2,b2.conjugate())]])
  sig_y=SIGMA_Y_GEN(4)
  rho2=sig_y*rho.conjugate()*sig_y
  E,W=eig(sqrt(sqrt(rho)*rho2*sqrt(rho)))
  E=sort(E)  
  return max([0.,E[3]-E[2]-E[1]-E[0]])
  

def H(self) :
  return self.conjugate().T

def GUE(dim) :
  A=randn(dim,dim) + 1j*randn(dim,dim)
  return sqrt(1/(2.*dim))*( A + H(A)) 
  
def RandState(n) :
  a=rand(n) + 1j*rand(n)
  return a/norm(a)

def BasisState(n) :
  b=1/sqrt(2)*(array([1,1]))
  a=b
  while n!=2 :
    a=kron(a,b)
    n=n/2
  return a

def BellState() :
  b=1/sqrt(2)*(array([1,0,0,1]))
  return b

def SIGMA_X_SUM(nc) :
  H=SIGMA_X
  i=nc
  Ic=2
  while i!=2 :
    H=kron(H,eye(2))+kron(eye(Ic),SIGMA_X)
    Ic+=2
    i=i/2
  return H

def SIGMA_Y_GEN(nc) :
  H=SIGMA_Y
  i=nc
  while i!=2 :
    H=kron(H,SIGMA_Y)
    i=i/2
  return H
  

def Hamiltonian_CBA_gen(dim_c,dim_b,dim_a,theta,gamma,lamb) :
  Ha=GUE(dim_a)
  Hb=GUE(dim_b)
  Hc=eye(dim_c)
  Vab=GUE(dim_a*dim_b)
  Vbc=GUE(dim_b*dim_c)
  Ia=eye(dim_a)
  Ib=eye(dim_b)
  Ic=eye(dim_c)
  
  Henv=cos(theta)*kron(Ha,Ib)+kron(Hb,Ia)+sin(theta)*gamma*Vab
  H=kron(Hc,eye(dim_a*dim_b))+ kron(Henv,Ic)+ lamb*kron(Vbc,Ia)
  return H


def Hamiltonian_CBA(dim_c,dim_b,dim_a,lamb) :
  Hc=SIGMA_X_SUM(dim_c)
  Hab=GUE(dim_a*dim_b)
  Vcb=GUE(dim_c*dim_b)
  Ia=eye(dim_a)
  Ib=eye(dim_b)
  Ic=eye(dim_c)
  
  H=kron(Hc,eye(dim_a*dim_b))+kron(Ic,Hab)+lamb*kron(Vcb,Ia)
  return H

def Hamiltonian_CBA_2special(dim_c,dim_b,dim_a,lamb) :
  Hc=SIGMA_X_SUM(dim_c)
  Hab=GUE(dim_a*dim_b)
  Vcb=GUE(2*dim_b)
  Ia=eye(dim_a)
  Ib=eye(dim_b)
  Ic=eye(dim_c)
  
  H=kron(Hc,eye(dim_a*dim_b))+kron(Ic,Hab)+lamb*kron(kron(Vcb,eye(2)),Ia)
  return H

