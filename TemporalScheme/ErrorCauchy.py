from numpy import array, zeros, linspace
from Esquemas import Euler, RK2, RK4, EulerI, CrankNicholson, Cauchy
from Conv import refine_mesh

############################################################################
############################### ERROR ######################################
############################################################################

def ErrorCauchy(F, U0, t, N, Esquema): # (Funcion, U0, t, N, Esquema)
    
    U1 = Cauchy(F, t, N, U0, Esquema)
    U2 = Cauchy(F, t, 2*N, U0, Esquema)
    
    Error = array(zeros((N, len(U0))))

    for i in range(0,N):
        Error[i,:] = (U2[2*i,:]-U1[i,:]) # E = (U2(dt/2)-U1(dt))

    return Error 


#def ErrorCauchy2(F, U0, t, N, Esquema): # (Funcion, U, U0, t, N, Esquema)
#    
#    U1 = Cauchy(F, t, N, U0, Esquema)
#    U2 = Cauchy(F, t, 2*N, U0, Esquema)
#    
#    Error = array(zeros((N, len(U0))))
#    if Esquema == Euler:
#        q = 1
#    elif Esquema == RK2:
#        q = 2
#    elif Esquema == RK4:
#        q = 4
#    elif Esquema == CrankNicholson: 
#        q = 2
#
#    for i in range(0,N):
#        Error[i,:] = (U2[2*i,:]-U1[i,:])/(1-1/(2**q)) # E = (U2(dt/2)-U1(dt))/(1-1/2**q)
#
#    return Error 