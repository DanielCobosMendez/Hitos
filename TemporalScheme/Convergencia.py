from numpy import array, zeros, linspace, log10
from TemporalScheme.ErrorCauchy import ErrorCauchy
from numpy.linalg import norm

def Conv(F, U0, t, N, Esquema):
    E = ErrorCauchy(F, U0, t, N, Esquema)
    q = array(zeros((N, len(U0))))
    logN = array(zeros((N, 1)))
    for i in range(4, N):
        q[i,:] = (log10(E[i,:])-log10(E[i-1,:]))/(log10(i)-log10(i-1))
        logN[i,0] = log10(i)
        print(q[i,:])
    
    return q, logN