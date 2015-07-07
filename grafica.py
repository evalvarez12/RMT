from pylab import *

lines=['o-','<-','p-','v-','D-','^-','h-','D-','>-','H-','d-','s-','x-','*-']
li=0


#for lamb in [.001,.01,.05,.1,.5] :
  #PS=load("g1_p_lambda="+str(lamb)+"_2q.npy")
  #plot(PS,lines[li],markevery=50,label="$\lambda="+str(lamb)+"$")
  #li+=1
  
  
  
#ylabel("$C$",fontsize=28)
#xlabel("$t$",fontsize=28)
#legend(loc='lower left',fontsize=25)

#yticks((1,0.9,0.8,0.7,0.6,0.5),('$1$','$0.9$','$0.8$','$0.7$','$0.6$','$0.5$'),fontsize=25)
#xticks((0,500,1000,1500,2000,2501),('$0$','$500$','$1000$','$1500$','$2000$','$2500$'),fontsize=25)

#axis([0,2500,.5,1])

#ax = axes([.55, .65, .33, .22])

#li=0
#for lamb in [.001,.01,.05,.1,.5] :
  #ains=load("g3_p_lambda="+str(lamb)+"_2q.npy")
  #plot(ains,lines[li],markevery=50)
  #li+=1
  ##xticks((0,pi/8,pi/4),('$0$','$\pi/8$','$\pi/4$'),fontsize=25)
  ##yticks((0.9,0.8,0.7,0.6,0.5),('$0.9$','$0.8$','$0.7$','$0.6$','$0.5$'),fontsize=25)
  #setp(ax)


#ylabel("$P$",fontsize=28)
#xlabel("$t$",fontsize=28)

#yticks((1,0.8,0.6,0.4,0.2),('$1$','$0.8$','$0.6$','$0.4$','$0.2$'),fontsize=25)
#xticks((0,500,1000,1500,2000,2501),('$0$','$500$','$1000$','$1500$','$2000$','$2500$'),fontsize=25)




#############  FIGURA   4   ##############

a=load("g4_p_Asize=2.npy")
b=load("g1_p_lambda=0.01_2q.npy")
c=load("g4_p_Asize=4.npy")
d=load("g4_p_Asize=5.npy")

plot(a,lines[0],markevery=50,label="$N_e=2^2$")
plot(b,lines[1],markevery=50,label="$N_e=2^3$")
plot(c,lines[2],markevery=50,label="$N_e=2^4$")
plot(d,lines[3],markevery=50,label="$N_e=2^5$")

yticks((1,0.9,0.8,0.7,0.6,0.5),('$1$','$0.9$','$0.8$','$0.7$','$0.6$','$0.5$'),fontsize=25)
xticks((0,500,1000,1500,2000,2501),('$0$','$500$','$1000$','$1500$','$2000$','$2500$'),fontsize=25)

axis([0,2500,.75,1])

ylabel("$P$",fontsize=28)
xlabel("$t$",fontsize=28)

legend(loc='upper right',fontsize=25)



##########  FIGURA  2   ################

#a=load("Pxlambda.npy")
#lamb = arange(0.001,1.005,0.005)
#plot(lamb,a,'k')

#yticks((1,0.9,0.8,0.7,0.6,0.5),('$1$','$0.9$','$0.8$','$0.7$','$0.6$','$0.5$'),fontsize=25)
#xticks((0.001,0.01,.1,1),('$10^{-3}$','$10^{-2}$','$10^{-1}$','$10^{0}$'),fontsize=25)

#axis([0.001,1,.5,1])
#xscale('log')
#ylabel("$P$",fontsize=28)
#xlabel("$\lambda$",fontsize=28)



show()
