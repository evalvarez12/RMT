from pylab import * 
import rmt


def H(self) :
  return self.conjugate().T

def GUE(dim) :
  A=randn(dim,dim) + 1j*randn(dim,dim)
  return sqrt(1/(4.*dim))*( A + H(A)) 
  
  



def Hamiltonian_CBA(dim_c,dim_b,dim_a,theta,gamma,lamb) :
  Ha=GUE(dim_a)
  Hb=GUE(dim_b)
  Hc=identity(dim_c)
  Vab=GUE(dim_a*dim_b)
  Vbc=GUE(dim_b*dim_c)
  Ia=identity(dim_a)
  Ib=identity(dim_b)
  Ic=identity(dim_c)
  
  Henv=cos(theta)*kron(Ha,Ib)+kron(Hb,Ia)+sin(theta)*gamma*Vab
  H=kron(Hc,identity(dim_a*dim_b))+ kron(Henv,Ic)+ lamb*kron(Vbc,Ia)
  return H


nc=2
na=8
nb=6
lamb=0.1

H=Hamiltonian_CBA(nc,nb,na,0,0,lamb)