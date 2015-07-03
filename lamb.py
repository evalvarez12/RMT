from pylab import * 
import rmt_py as rmt




def one_run_fixed_t(Ham_func,nc,nb,na,lamb,time) :
  H=Ham_func(nc,nb,na,lamb)
  E,W=eig(H)
  state=kron(rmt.BasisState(nc),kron(rmt.RandState(nb),rmt.RandState(na)))
  stateW=dot(W.conjugate().transpose(),state)
  purities=[]
  state_t=dot(W,exp(-1j*E*time)*stateW)
  return rmt.purity_C(state_t).real



nc=2
na=2**10
nb=2**6


time=1000

lines=['o-','<-','p-','v-','+-','^-','h-','D-','>-','H-','d-','s-','x-','*-']
li=0
ave=50
Ps=[]
for lamb in arange(0.001,1.005,0.005) :
  p_ave=0
  for j in range(ave) :
    p=one_run_fixed_t(rmt.Hamiltonian_CBA,nc,nb,na,lamb,time)
    p_ave+=p
    
  print lamb  
  p_ave=p_ave/double(ave)
  Ps.append(p_ave)


save("Pxlambda",Ps)
plot(arange(0.001,1.005,0.005),Ps,markevery=20)

  
  
ylabel("$P$",fontsize=28)
xlabel("$\lambda$",fontsize=28)

show()

