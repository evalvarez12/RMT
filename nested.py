from pylab import * 


SIGMA_X=array([[0,1],[1,0]])
SIGMA_Y=array([[0,-1j],[1j,0]])
SIGMA_Z=array([[1,0],[0,-1]])

def purity_C(state) :
  l=size(state)
  a=state[:l/2]
  b=state[l/2:]
  rho=array([[dot(a,a.conjugate()),dot(a,b.conjugate())],[dot(b,a.conjugate()),dot(b,b.conjugate())]])
  return trace(dot(rho,rho))

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
  Hc=SIGMA_X
  Hab=GUE(dim_a*dim_b)
  Vcb=GUE(dim_c*dim_b)
  Ia=eye(dim_a)
  Ib=eye(dim_b)
  Ic=eye(dim_c)
  
  H=kron(Hc,eye(dim_a*dim_b))+kron(Ic,Hab)+lamb*kron(Vcb,Ia)
  return H

nc=2
na=10
nb=6
lamb=.1

time=1000

H=Hamiltonian_CBA(nc,nb,na,lamb)
E,W=eig(H)

state=kron(BasisState(nc),kron(RandState(nb),RandState(na)))

stateW=dot(W.conjugate().transpose(),state)

purs=[]

for i in range(time) :
  state_t=dot(W,exp(-1j*E*i)*stateW)
  print real(purity_C(state_t))
  purs+=[real(purity_C(state_t))]
  
plot(purs)  
#show()

