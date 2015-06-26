#RMT in julia
#Eduardo Villaseñor - evalvarez12@gmail.com

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
  while n!=2 inv
    a=kron(a,b)
    n=n/2
  end
  return a
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


