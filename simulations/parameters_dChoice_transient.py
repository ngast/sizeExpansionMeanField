import numpy as np

def fixedPoint(n,rho,d):
    return( np.array([rho**( (d**(i+1)-1)/(d-1)) for i in range(n)],dtype=np.float32) )

def defineDrift(n,rho,d,evaluate_at=None):
    def computeF(x):
        F = np.zeros((n),dtype=np.float32)
        for i in range(1,n-1):
            F[i] = rho*(x[i-1]**d-x[i]**d)+x[i+1]-x[i]
        F[0] = rho*(1-x[0]**d)+x[1]-x[0]
        F[n-1] = rho*(x[n-2]**d-x[n-1]**d)-x[n-1]
        return(F)
    if evaluate_at is None:
        return(computeF)
    else:
        return(computeF(evaluate_at))

def defineDriftDerivativeQ(n,rho,d,evaluate_at=None):
    def computeFp(x):
        M = np.zeros((n,n), dtype=np.float32)
        for i in range(n):
            M[i,i] = -1 - d*rho*x[i]**(d-1)
        for i in range(n-1):
            M[i,i+1] = 1;
            M[i+1,i] = d*rho*x[i]**(d-1)
        return(M)
    
    def computeFpp(x):
        fPP = np.zeros((n,n,n), dtype=np.float32)
        for i in range(n):
            fPP[i][i][i] = -d*(d-1)*rho*x[i]**(d-2)
        for i in range(n-1):
            fPP[i+1][i][i] = d*(d-1)*rho*x[i]**(d-2)
        return(fPP)
    
    def computeQ(x):
        Q = np.zeros((n,n),dtype=np.float32)
        for i in range(1,n-1):
            Q[i,i] = rho*(x[i-1]**d-x[i]**d) + x[i]-x[i+1]
        Q[0,0] = rho*(1-x[0]**d) + x[0]-x[1]
        Q[n-1,n-1] = rho*(x[n-2]**d-x[n-1]**d) + x[n-1]
        return(Q)
    if evaluate_at is None:
        return(computeFp,computeFpp,computeQ)
    else:
        x = evaluate_at
        return(computeFp(x),computeFpp(x),computeQ(x))

def defineDriftSecondDerivativeQderivativesR(n,rho,d,evaluate_at=None):        
    def computeFppp(x):
        fPPP=np.zeros((n,n,n,n),dtype=np.float32)
        for i in range(n):
            fPPP[i][i][i][i] = -d*(d-1)*(d-2)*rho*x[i]**(d-3) if d >= 3 else 0
        for i in range(n-1):
            fPPP[i+1][i][i][i] = d*(d-1)*(d-2)*rho*x[i]**(d-3) if d >= 3 else 0
        return(fPPP)

    def computeFpppp(x):
        fPPPP=np.zeros((n,n,n,n,n),dtype=np.float32)
        for i in range(n):
            fPPPP[i][i][i][i][i] = -d*(d-1)*(d-2)*(d-3)*rho*x[i]**(d-4) if d >= 4 else 0
        for i in range(n-1):
            fPPPP[i+1][i][i][i][i] = d*(d-1)*(d-2)*(d-3)*rho*x[i]**(d-4) if d >= 4 else 0
        return(fPPPP)
    def computeQp(x):
        Qp = np.zeros((n,n,n),dtype=np.float32)
        for i in range(n):
            Qp[i,i,i] = 1-d*rho*x[i]**(d-1)
            if i>0: Qp[i,i,i-1] = d*rho*x[i-1]**(d-1)
            if i<n-1: Qp[i,i,i+1] = -1
        return(Qp)
    def computeQpp(x):
        Qpp = np.zeros((n,n,n,n),dtype=np.float32)
        for i in range(n):
            Qpp[i,i,i,i] = -d*(d-1)*rho*x[i]**(d-2)
            if i>0: Qpp[i,i,i-1,i-1] = d*(d-1)*rho*x[i-1]**(d-2)
        return(Qpp)
    computeF = defineDrift(n,rho,d,evaluate_at=None)
    def computeR(x):
        F = computeF(x)
        R = np.zeros((n,n,n),dtype=np.float32)
        for i in range(n):
            R[i,i,i] = F[i]
        return(R)
    if evaluate_at is None:
        return(computeFppp,computeFpppp,computeQp,computeQpp,computeR)
    else:
        x = evaluate_at
        return(computeFppp(x),computeFpppp(x),computeQp(x),computeQpp(x),computeR(x))
