from numpy import array, zeros, linspace, log10, polyfit
from TemporalScheme.ErrorCauchy import ErrorCauchy
from numpy.linalg import norm

def Conv(F, U0, t, N0, Esquema, T):
    """
    Inputs:
        F : Función del problema de Cauchy
        U0 : Vector de estado inicial
        t : Tiempo total
        N0 : Número de particiones iniciales
        Esquema: Esquema numérico

    """

    np = 6  # Número de puntos de la regresión, si se sube más tarda MUCHO en converger
    logE = zeros(np)
    logN = zeros(np)
    print(t)
    N = len(t) - 1
    t1 = t
    
    for i in range(np):
        E = ErrorCauchy(F, U0, T, N, Esquema)  # Asumiendo que Error devuelve U_1 y Error
        logE[i] = log10(norm(E[N, :])) 
        logN[i] = log10(float(N))
        N = 2*N
        t1 = linspace(t[0], t[-1], N+1) 

    y = logE[ logE > -12 ]
    x = logN[ 0:len(y) ]
    Order, b = polyfit(x, y, 1) # Regresión lineal para encontrar la pendiente de la recta que mejor se ajusta a los datos

    #ptosgraf = 8
    #logE = zeros(ptosgraf)
    #logN = zeros(ptosgraf)
    #N = N0
    #for i in range(0,ptosgraf):
    #    
    #    E = ErrorCauchy(F, U0, t, N, Esquema)
    #    logE[i] = log10(norm(E[-1,:]))
    #    logN[i] = log10(N)
    #    N = N0**i

    return logE, logN, Order

#def Convergence(U0, F, Error, Problema, Esquema, t):  
#    np = 6  # Número de puntos de la regresión, si se sube más tarda MUCHO en converger
#    logE = zeros(np)
#    logN = zeros(np)
#    N = len(t) - 1
#    t1 = t
#    
#    for i in range(np):
#        U, Er = Error(U0, F, Problema, Esquema, t1)  # Asumiendo que Error devuelve U_1 y Error
#        logE[i] = log10(norm(Er[N, :])) 
#        logN[i] = log10(float(N))
#        N = 2*N
#        t1 = linspace(t[0], t[-1], N+1)  
#        
#    y = logE[ logE > -12 ]
#    x = logN[ 0:len(y) ]
#    Order, b = polyfit(x, y, 1) # Regresión lineal para encontrar la pendiente de la recta que mejor se ajusta a los datos
#    # print("Order =", Order, "b =", b)
#    return logN, logE, Order