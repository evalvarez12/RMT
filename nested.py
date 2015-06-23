from pylab import * 
import rmt

def Hamiltonian_CBA(dim_c,dim_b,dim_a,theta,gamma,lamb) :
  Ha=GUE(dim_a)
  Hb=GUE(dim_b)
  Hc=identity(dim_c)
  Vab=GUE(dim_a*dim_b)
  Vbc=GUE(dim_b*dim_a)
  Ia=identity(dim_a)
  Ib=identity(dim_b)
  Ic=identity(dim_c)
  
  Henv=cos(theta)*kron(Ha,Ib)+kron(Hb,na)+sin(theta)*gamma*Vab
  H=kron(Hc,identity(dim_a*dim_b))+ kron(Henv,Ic)+ lamb*kron(Vbc,Ia)
  return H


nc=1
na=8
nb=6
theta=0

H=Hamiltonian_CBA(nc,nb,na,)