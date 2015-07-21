# main.jl

@everywhere include("rmt.jl")


function one_run(Hamiltonian_func,nb,na,lamb,time) 
  ave=50
  P_ave=zeros(time)
  for i in range(0,ave)
    state=kron(BasisState(2),kron(RandState(nb),RandState(na)))
    H=Hamiltonian_func(nb,na,lamb)
    E,W=eig(H)
    state=(W')*state
    purs=[]
    for i in range(0,time)
      state_t=W*(exp(-1*im*E*i).* state)
      append!(purs,[real(purity(state_t))])
    end
    P_ave+=purs
  end  
  return P_ave/Float64(ave)
end


na=2^8
nb=2^4
lamb=.01

time=2000

p=one_run(Hamiltonian_CBA,nb,na,lamb,time) 
writedlm("g1_p_lamb=.01_na=2^8_nb=2^4.dat",p)


# plot(purs)  