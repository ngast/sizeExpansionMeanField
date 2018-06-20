import parameters_2choice_transient as twoChoice
import rmf_tool.src.rmf_tool as rmf
import numpy as np
import matplotlib.pylab as plt 

class twoChoiceRefined(rmf.DDPP):
    def __init__(self,n,rho):
        super(twoChoiceRefined,self).__init__()
        self.add_transition(np.zeros(n),lambda x:x)
        self.n = n
        self.rho=rho
        self.set_initial_state(np.zeros(n))

    def fixed_point(self):
        return(twoChoice.fixedPoint(self.n,self.rho,2))
        
    def defineDrift(self,evaluate_at=None):
        return(twoChoice.defineDrift(self.n,self.rho,evaluate_at))

    def defineDriftDerivativeQ(self,evaluate_at=None):
        return(twoChoice.defineDriftDerivativeQ(self.n,self.rho,evaluate_at))

    def defineDriftSecondDerivativeQderivativesR(self,evaluate_at=None):
        return(twoChoice.defineDriftSecondDerivativeQderivativesR(self.n,self.rho,evaluate_at))

K=20
rho=0.9

myModel = twoChoiceRefined(K,rho)
X,V,_ = myModel.meanFieldExpansionSteadyState(order=1)
X,V,A,_ = myModel.meanFieldExpansionSteadyState(order=2)
for N in [10,20,30]:
    print(np.sum(X),np.sum(X+V/N), np.sum(X+V/N+A/N**2))

x0 = np.zeros(K)
x0[0:3] = [1,1,.8]
myModel.set_initial_state(x0)

N=10
T,X,V,_ = myModel.meanFieldExpansionTransient(order=1,time=100)
plt.clf()
plt.plot(T,np.sum(X,1))
plt.plot(T,np.sum(X+V/N,1),'--')

import rmf_tool.misc.jsqD_simulate.average_valueJSQ as jsqSimu
Tsimu,Xsimu = jsqSimu.loadTransientSimu(N,2,rho,2.8)
plt.plot(Tsimu,np.sum(Xsimu,1)-1)
