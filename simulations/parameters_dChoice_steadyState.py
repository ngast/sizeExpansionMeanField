import numpy as np

def fixedPoint(parameters):
    (n,d,rho)=parameters
    return( np.array([rho**( (d**(i+1)-1)/(d-1)) for i in range(n)],dtype=np.float32) )
def computeFp(parameters):
    (n,d,rho)=parameters
    pi = fixedPoint(parameters)
    M = np.zeros((n,n), dtype=np.float32)
    for i in range(n):
        M[i,i] = -1 - d*rho*(pi[i])**(d-1)
    for i in range(n-1):
        M[i,i+1] = 1;
        M[i+1,i] = d*rho*(pi[i])**(d-1)
    return(M)
def computeFpp(parameters):
    (n,d,rho)=parameters
    pi = fixedPoint(parameters)
    fPP = np.zeros((n,n,n), dtype=np.float32)
    for i in range(n):
        fPP[i][i][i] = -d*(d-1)*rho*(pi[i])**(d-2)
        if i<n-1:
            fPP[i+1][i][i] = d*(d-1)*rho*(pi[i])**(d-2)
    return(fPP)
def computeFppp(parameters):
    (n,d,rho)=parameters
    if d<=2:
        return(np.zeros((n,n,n,n),dtype=np.float32))
    else:
        pi = fixedPoint(parameters)
        Fppp = np.zeros((n,n,n,n),dtype=np.float32)
        for i in range(n):
            Fppp[i,i,i,i] = -d*(d-1)*(d-2)*rho*(pi[i])**(d-3)
            if i<n-1:
                Fppp[i+1,i,i,i] = d*(d-1)*(d-2)*rho*(pi[i])**(d-3)
        return(Fppp)
def computeFpppp(parameters):
    (n,d,rho)=parameters
    if d<=3:
        return(np.zeros((n,n,n,n,n),dtype=np.float32))
    else:
        pi = fixedPoint(parameters)
        Fpppp = np.zeros((n,n,n,n,n),dtype=np.float32)
        for i in range(n):
            Fpppp[i,i,i,i,i] = -d*(d-1)*(d-2)*(d-3)*rho*(pi[i])**(d-4)
            if i<n-1:
                Fpppp[i+1,i,i,i,i] = d*(d-1)*(d-2)*(d-3)*rho*(pi[i])**(d-4)
        return(Fpppp)

def computeQ(parameters):
    (n,d,rho)=parameters
    pi = fixedPoint(parameters)
    Q = np.zeros((n,n),dtype=np.float32)
    Q[0,0] = rho*(1 - pi[0]**d)+(pi[0]-pi[1])
    for i in range(1,n-1):
        Q[i,i] = rho*(pi[i-1]**d - pi[i]**d)+(pi[i]-pi[i+1])
    Q[n-1,n-1] = rho*(pi[n-2]**d - pi[n-1]**d)+(pi[n-1])
    return(Q)

def computeR(parameters):
    (n,d,rho)=parameters
    R = np.zeros((n,n,n),dtype=np.float32)
    return(R)

def computeQp(parameters):
    (n,d,rho)=parameters
    pi = fixedPoint(parameters)
    Qp = np.zeros((n,n,n),dtype=np.float32)
    for i in range(n):
        Qp[i,i,i] = 1-d*rho*(pi[i])**(d-1)
        if i<n-1:
            Qp[i,i,i+1] = -1
            Qp[i+1,i+1,i] = d*rho*(pi[i])**(d-1)
    return(Qp)

def computeQpp(parameters):
    (n,d,rho)=parameters
    pi = fixedPoint(parameters)    
    Qpp = np.zeros((n,n,n,n),dtype=np.float32)
    for i in range(n):
        Qpp[i,i,i,i] = -d*(d-1)*rho*(pi[i])**(d-2)
        if i<n-1:
            Qpp[i+1,i+1,i,i] = d*(d-1)*rho*(pi[i])**(d-2)
    return(Qpp)
