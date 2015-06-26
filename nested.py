from pylab import * 
import rmt_py as rmt



#def one_run(Ham_func,nc,nb,nc,lamb,time) :
  #H=Ham_func(nc,nb,na,lamb)
  #E,W=eig(H)
  #state=kron(rmt.BasisState(nc),kron(rmt.RandState(nb),rmt.RandState(na)))
  #stateW=dot(W.conjugate().transpose(),state)
  #for i in range(time) :
    #state_t=dot(W,exp(-1j*E*i)*stateW)
    ##print real(rmt.purity_C(state_t))
    #purs.append(real(rmt.purity_C(state_t)))

nc=2
na=10
nb=6
lamb=.1

time=1000

PS=zeros(time)

ave=50
for lamb in [.001,.01,.05,.1,.5] :
  for j in range(ave) :
    purs=[]
    H=rmt.Hamiltonian_CBA(nc,nb,na,lamb)
    E,W=eig(H)
    state=kron(rmt.BasisState(nc),kron(rmt.RandState(nb),rmt.RandState(na)))
    stateW=dot(W.conjugate().transpose(),state)
    for i in range(time) :
      state_t=dot(W,exp(-1j*E*i)*stateW)
      #print real(rmt.purity_C(state_t))
      purs.append(real(rmt.purity_C(state_t)))
    PS+=purs
    print j
    
  PS=PS/double(ave)
  plot(PS,label="$"+str(lamb)+"$")


label
legend()
show()

