from pylab import * 
import rmt_py as rmt



def one_run(Ham_func,nc,nb,na,lamb,time) :
  H=Ham_func(nc,nb,na,lamb)
  E,W=eig(H)
  state=kron(rmt.BasisState(nc),kron(rmt.RandState(nb),rmt.RandState(na)))
  #state=kron(rmt.BellState(),kron(rmt.RandState(nb),rmt.RandState(na)))
  stateW=dot(W.conjugate().transpose(),state)
  purities=[]
  #concurrences=[]
  for i in range(time) :
    state_t=dot(W,exp(-1j*E*i)*stateW)
    #print real(rmt.purity_C(state_t))
    purities.append(rmt.purity_C(state_t).real)
    #concurrences.append(rmt.concurence(state_t))
  #return [array(purities),array(concurrences)]
  return array(purities)




nc=2
na=2**6
nb=2**3


time=2500

lines=['o-','<-','p-','v-','+-','^-','h-','D-','>-','H-','d-','s-','x-','*-']
li=0
ave=50
for lamb in [.001,.01,.05,.1,.5] :
  PS=zeros(time)
  #CS=zeros(time)
  for j in range(ave) :
    purs=one_run(rmt.Hamiltonian_CBA,nc,nb,na,lamb,time)
    PS+=purs
    #CS+=concs
    
  print lamb  
  PS=PS/double(ave)
  #CS=CS/double(ave)
  save("g1_p_lambda="+str(lamb)+"_2q",PS)
  #save("g3_c_lambda="+str(lamb)+"_2q",CS)
  #plot(PS,lines[li],markevery=50,label="$\lambda="+str(lamb)+"$")
  #li+=1
  
  
#ylabel("$P$",fontsize=28)
#xlabel("$t$",fontsize=28)
#legend(loc='lower left',fontsize=25)
#show()

