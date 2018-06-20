from rmf_tool.src.refinedRefined_fixedPoint import *
from parameters_dChoice_steadyState import *
import rmf_tool.misc.jsqD_simulate.average_valueJSQ as jsqSimu

myN=[10,20,30,50,100]
for N in myN:
    print(' ',N,'   (simu) ',end='\t| ')
print()

from numpy import set_printoptions
set_printoptions(precision=4)

K=20 # queue lenght are bounded by 20

for d in [2,3,4]:
    print('d=',d)
    for rho in [0.7,0.9,0.95]:
        parameters = (K,d,rho)
        pi=fixedPoint(parameters)
        Fp=computeFp(parameters)
        Fpp=computeFpp(parameters)
        Fppp=computeFppp(parameters)
        Fpppp=computeFpppp(parameters)
        Q=computeQ(parameters)
        Qp=computeQp(parameters)
        Qpp=computeQpp(parameters)
        R=computeR(parameters)
        pi,V,A, VWABCD = computePiVA(pi,Fp,Fpp,Fppp,Fpppp,Q,Qp,Qpp,R)
        for N in myN:
            simu = jsqSimu.loadSteadyStateAverageQueueLength(N,d,rho)
            print('{0:.4f} ({1:.4f})'.format(sum(pi+V/N+A/N**2),simu),end='\t| ')
        print('')

import matplotlib.pyplot as plt

myRho = np.array([0.7,0.75,0.8,0.85,0.9,0.95,0.97,0.99])
myD = [2,3,4,5]
myPi = np.zeros((len(myRho),4))
myV = np.zeros((len(myRho),4))
for iD,d in enumerate(myD):
    for i,rho in enumerate(myRho):
        parameters = (K,d,rho)
        pi=fixedPoint(parameters)
        Fp=computeFp(parameters)
        Fpp=computeFpp(parameters)
        #Fppp=computeFppp(parameters)
        #Fpppp=computeFpppp(parameters)
        Q=computeQ(parameters)
        #Qp=computeQp(parameters)
        #Qpp=computeQpp(parameters)
        #R=computeR(parameters)
        pi,V,_ = computePiV(pi,Fp,Fpp,Q)
        myPi[i,iD] = np.sum(pi)
        myV[i,iD] = np.sum(V)

plt.figure()
for iD,d in enumerate(myD):
    plt.plot(myRho, myV[:,iD]*(1-myRho)/myRho**2 )
plt.legend(['d={}'.format(d) for d in myD])
plt.show()
