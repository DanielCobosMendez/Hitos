from numpy import array, zeros, linspace
from TemporalScheme.Esquemas import Euler, RK2, RK4, EulerI, CrankNicholson

############################################################################
############################### ERROR ######################################
############################################################################

def Error(U1, U2, U0, N, Esquema):
    Error = array(zeros((N, len(U0))))
    if Esquema == Euler:
        q = 1
    elif Esquema == RK2:
        q = 2
    elif Esquema == RK4:
        q = 4
    elif Esquema == CrankNicholson: 
        q = 2

    for i in range(0,N):
        Error[i,:] = (U2[2*i,:]-U1[i,:])/(1-1/(2**q)) # E = (U2(dt/2)-U1(dt))/(1-1/2**q)

    return Error