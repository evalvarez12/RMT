# main.jl

include("rmt.jl")




nc=2
na=10
nb=6
lamb=.1

time=10

H=Hamiltonian_CBA(nc,nb,na,lamb)
E,W=eig(H)

state=kron(BasisState(nc),kron(RandState(nb),RandState(na)))

stateW=(W')*state

purs=[]

l=size(state)[1]
l=round(Int,l/2)

for i in range(0,time) 
  state_t=W*(exp(-1*im*E*i).* stateW)
  println(real(purity(state_t)))
  #CALCULA 
#   a=state_t[1:l]
#   b=state_t[l+1:end]
#   rho=[dot(a,conj(a)) dot(a,conj(b)) ; dot(b,conj(a)) dot(b,conj(b))]
#   p=trace(rho*rho)
#   println(real(p))
  
  #append!(purs,[real(purity_C(state_t))])
end
  
# plot(purs)  