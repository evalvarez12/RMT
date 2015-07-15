# main.jl

@everywhere include("rmt.jl")


function one_run(Hamiltonian_func,nc,nb,na,lamb,time) 
  state=kron(BasisState(nc),kron(RandState(nb),RandState(na)))
  H=Hamiltonian_func(nc,nb,na,lamb)
  E,W=eig(H)
  state=(W')*state
  purs=[]
  for i in range(0,time)
    state_t=W*(exp(-1*im*E*i).* state)
    append!(purs,[real(purity(state_t))])
  end
  return purs
end

nc=2
na=2^8
nb=2^4
lamb=.1

time=1

p=one_run(Hamiltonian_CBA,nc,nb,na,lamb,time) 
writedlm("pursjl.dat",p)


# plot(purs)  