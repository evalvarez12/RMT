#RMT in julia
#Eduardo Villaseñor - evalvarez12@gmail.com

#kron(a,b) a queda al final del total

SIGMA_X=[0 1 ; 1 0]
SIGMA_Y=[ 0 -im*1 ; im*1 0]
SIGMA_Z=[1 0 ; 0 -1]


function GUE(n::Int) 
  A = randn(n, n) + im*randn(n, n)
  normalization = sqrt(2*n)
  return (A + A') / normalization
end
  
 
function RandState(n::Int)
  a=rand(n)*2-1 + im*rand(n)*2-1
  return a/norm(a)
end

function BasisState(n::Int)
  b=1/sqrt(2)*[1,1]
  a=b
  while n!=2 
    a=kron(a,b)
    n=n/2
  end
  return a
end

function BellState()
  b=1/sqrt(2)*[1,0,0,1]
  return b
end

function SIGMA_X_SUM(nc) 
  H=SIGMA_X
  i=nc
  Ic=2
  while i!=2 
    H=kron(H,eye(2))+kron(eye(Ic),SIGMA_X)
    Ic+=2
    i=i/2
  end
  return H
end

function SIGMA_Y_GEN(nc) 
  H=SIGMA_Y
  i=nc
  while i!=2 
    H=kron(H,SIGMA_Y)
    i=i/2
  end  
  return H
end  

function Hamiltonian_CBA(nc::Int,nb::Int,na::Int,lamb::Float64)
  Hc=SIGMA_X
  Hab=GUE(na*nb)
  Vcb=GUE(nc*nb)
  Ia=eye(na)
  Ib=eye(nb)
  Ic=eye(nc)
  
  H=kron(Hc,eye(na*nb))+kron(Ic,Hab)+lamb*kron(Vcb,Ia)
  return H
end

function Hamiltonian_CBA_2special(nc::Int,nb::Int,na::Int,lamb::Float64)
  Hc=SIGMA_X_SUM(nc)
  Hab=GUE(na*nb)
  Vcb=GUE(2*nb)
  Ia=eye(na)
  Ib=eye(nb)
  Ic=eye(nc)
  
  H=kron(Hc,eye(na*nb))+kron(Ic,Hab)+lamb*kron(kron(Vcb,eye(2)),Ia)
  return H
end


function purity(state) 
  l=size(state)[1]
  l=round(Int,l/2)
  a=state[1:l]
  b=state[l+1:end]
  rho=[dot(a,a) dot(b,a) ; dot(a,b) dot(b,b)]
  p=real(trace(rho*rho))
  return p
end

function purity_concurrence(state)
  l=size(state)[1]
  l=round(Int,l/2)
  l2=round(Int,l/2)
  a=state[1:l]
  a1=a[1:l2]
  a2=a[l2+1:end]
  b=state[l+1:end]
  b1=b[1:l2]
  b2=b[l2+1:end]
  rho=[dot(a1,a1) dot(a2,a1) dot(b1,a1) dot(b2,a1) ; dot(a1,a2) dot(a2,a2) dot(b1,a2) dot(b2,a2) ; dot(a1,b1) dot(a2,b1) dot(b1,b1) dot(b2,b1); dot(a1,b2) dot(a2,b2) dot(b1,b2) dot(b2,b2)]
  p=real(trace(rho*rho))
  sig_y=SIGMA_Y_GEN(4)
  rho2=sig_y*conj(rho)*sig_y
  E,W=eig(rho*rho2)
  E=sort(real(sqrt(E)))
  c=max(0.,E[4]-E[3]-E[2]-E[1])
  return p,c
end
