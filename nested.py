from pylab import * 
import rmt_py as rmt



def one_run(Ham_func,nc,nb,na,lamb,time) :
  H=Ham_func(nc,nb,na,lamb)
  E,W=eig(H)
  state=kron(rmt.BasisState(nc),kron(rmt.RandState(nb),rmt.RandState(na)))
  stateW=dot(W.conjugate().transpose(),state)
  purities=[]
  for i in range(time) :
    state_t=dot(W,exp(-1j*E*i)*stateW)
    #print real(rmt.purity_C(state_t))
    purities.append(rmt.purity_C(state_t).real)
  return array(purities)

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


time=2500

lines=['o-','<-','p-','v-','+-','^-','h-','D-','>-','H-','d-','s-','x-','*-']
li=0
ave=50
for lamb in [.001,.01,.05,.1,.5] :
  PS=zeros(time)
  for j in range(ave) :
    purs=one_run(rmt.Hamiltonian_CBA,nc,nb,nc,lamb,time)
    PS+=purs
    
  print lamb  
  PS=PS/double(ave)
  plot(PS,lines[ls],markevery=50,label="$\lambda="+str(lamb)+"$",)
  ls+=1
  
  
ylabel("$P$",fontsize=28)
xlabel("$t$",fontsize=28)
legend(loc='lower left'fontsize=25)
show()

