# main.jl

include("rmt.jl")

nc=2
na=10
nb=6
lamb=.1

time=1000

H=Hamiltonian_CBA(nc,nb,na,lamb)
E,W=eig(H)

state=kron(BasisState(nc),kron(RandState(nb),RandState(na)))

stateW=(W')*state

purs=[]

for i in range(1,time) :
  state_t=W*exp(-1*im*E*i*stateW)
#   print real(purity_C(state_t))
  append!(purs,[real(purity_C(state_t))])
end
  
# plot(purs)  