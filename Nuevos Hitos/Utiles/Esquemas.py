from scipy.optimize import newton
from numpy import array, zeros, linspace, log10, polyfit, sqrt
from numpy.linalg import norm

############################################################################
######################## MÉTODO DE EULER ###################################
############################################################################
"""
INPUTS:
    U: Vector de estado U(n)
    dt: Valor equiespaciado temporal
    F: Función del problema de Cauchy
    t: Paso temporal
OUTPUTS:
    U + dt * F(U, t): Vector de estado U(n+1) según el método de Euler
"""
def Euler(U, dt, F, t): # U(n+1) = U(n) + dt * F

    return U + dt * F(U, t)

############################################################################
############################### MÉTODO RK2 #################################
############################################################################
"""
INPUTS:
    U: Vector de estado U(n)
    dt: Valor equiespaciado temporal
    F: Función del problema de Cauchy
    t: Paso temporal
OUTPUTS:
    U + dt/2 * (K1_RK2 + K2_RK2): Vector de estado U(n+1) según el método de RK2
"""
def RK2(U, dt, F, t): # U(n+1) = U(n) + dt/2 * (K1 + K2)
    K1_RK2 = F(U, t)
    K2_RK2 = F(U + K1_RK2*dt, t + dt)
    return U + dt/2 * (K1_RK2 + K2_RK2)

############################################################################
############################### MÉTODO RK4 #################################
############################################################################
"""
INPUTS:
    U: Vector de estado U(n)
    dt: Valor equiespaciado temporal
    F: Función del problema de Cauchy
    t: Paso temporal
OUTPUTS:
    U + dt/6 * (K1_RK4 + 2 * K2_RK4 + 2 * K3_RK4 + K4_RK4): Vector de estado U(n+1) según el método de RK4
"""
def RK4(U, dt, F, t): # U(n+1) = U(n) + dt/6 * (K1 + 2*K2 + 2*K3 + K4)
    K1_RK4 = F(U, t)
    K2_RK4 = F(U + K1_RK4*dt/2, t + dt/2)
    K3_RK4 = F(U + K2_RK4*dt/2, t + dt/2)
    K4_RK4 = F(U + K3_RK4*dt, t + dt)
    return U + dt/6 * (K1_RK4 + 2 * K2_RK4 + 2 * K3_RK4 + K4_RK4)

############################################################################
######################## MÉTODO DE EULER INVERSO ###########################
############################################################################
"""
INPUTS:
    U: Vector de estado U(n)
    dt: Valor equiespaciado temporal
    F: Función del problema de Cauchy
    t: Paso temporal
OUTPUTS:
    newton(G, U, maxiter=500): Vector de estado U(n+1) según el método de Euler Inverso
"""
def EulerI(U, dt, F, t):
    def G(X):
        return X-U-dt*F(X, t)
    return newton(G, U, maxiter=1000)

############################################################################
###################### MÉTODO DE CRANK-NICOLSON ############################
############################################################################
"""
INPUTS:
    U: Vector de estado U(n)
    dt: Valor equiespaciado temporal
    F: Función del problema de Cauchy
    t: Paso temporal
OUTPUTS:
    newton(G, U, maxiter=500): Vector de estado U(n+1) según el método de Crank-Nicolson
"""
def CrankNicolson(U, dt, F, t):
    def G(X):
        return X-U-dt/2*(F(X, t)+F(U,t))
    return newton(G, U, maxiter=500)

############################################################################
######################## MÉTODO DE LEAP-FROG ###############################
############################################################################
"""
INPUTS:
    U: Vector de estado U(n-1)
    dt: Valor equiespaciado temporal
    F: Función del problema de Cauchy
    t: Paso temporal
    U1: Vector de estado U(n)
OUTPUTS:
    U(n+1) = U(n-1)+2*dt*F(U(n),t)
"""
def LeapFrog(U0, dt, F, t, U1): # U(n+1) = U(n-1) + 2 * dt * F
    return U0 + 2*dt * F(U1, t)

############################################################################
########################### REFINADO MALLA #################################
############################################################################
"""
INPUTS:
    t1: Malla temporal inicial
OUTPUTS:
    t2: Malla temporal refinada
"""
def refine_mesh(t1): 
      
      N = len(t1) - 1  
      t2 = zeros( 2*N +1) 

      for i in range(0,N): 
           t2[2*i] =  t1[i]
           t2[2*i+1] = ( t1[i]  + t1[i+1] )/ 2 
      t2[2*N] = t1[N]      

      return t2 

############################################################################
########################### PROBLEMA DE CAUCHY #############################
############################################################################
"""
INPUTS:
    F: Función del problema de Cauchy
    dt: Vector del equiespaciado temporal
    U0: Vector de estado inicial
    Esquema: Esquema temporal utilizado
    t: Periodo
    N: Particiones
OUTPUTS:
    U: Vector de estado según el esquema numérico y los parámetros escogidos
"""
def Cauchy(F, dt, U0, Esquema, T, N):
    U = array(zeros((len(dt), len(U0)))) # Definición de U
    U[0,:] = U0 # Asignación del vector de estado inicial
    if Esquema != LeapFrog:
        for i in range(len(dt)-1):
            U[i+1, :] = Esquema(U[i,:], dt[i+1]-dt[i], F, T) # Llamada al esquema numérico, e integración por cada paso temporal
    else:
        for i in range(1, len(dt)-1):
            U[1,:] = U[0,:]+dt[1]-dt[0]*F(U0, T)
            U[i+1, :] = Esquema(U[i,:], dt[i+1]-dt[i], F, T, U[i,:])

    return U

############################################################################
#################### PROBLEMA DE CAUCHY + ERROR ############################
############################################################################
"""
INPUTS:
    F: Función del problema de Cauchy
    dt: Vector del equiespaciado temporal
    U0: Vector de estado inicial
    Esquema: Esquema temporal utilizado
    T: Periodo
    N: Particiones
OUTPUTS:
    U1: Vector de estado según el esquema numérico y los parámetros escogidos
    Error: Error entre vectores de estado al cambiar el número de puntos de interpolación
"""
def ErrorCauchy(F, dt, U0, Esquema, T, N):
    
    t1 = dt
    t2 = refine_mesh(t1)
    N2 = len(t2)-1

    U1 = Cauchy(F, t1, U0, Esquema, T, N)
    U2 = Cauchy(F, t2, U0, Esquema, T, N2)

    Error = array(zeros((len(t1), len(U0))))

    for i in range(0, len(t1)):
        Error[i, :] = (U2[2*i, :] - U1[i, :])

    return Error, U1

############################################################################
########################## CONVERGENCIA ####################################
############################################################################
"""
INPUTS:
    F: Función del problema de Cauchy
    dt: Vector del equiespaciado temporal
    U0: Vector de estado inicial
    Esquema: Esquema temporal utilizado
    T: Periodo
    N: Particiones
OUTPUTS:
    y: valores de la convergencia
    x: valores de N para los valores de la convergencia
    order: orden del esquema numérico
"""
def Conv(F, dt, U0, Esquema, T, N):
    m = 6 # Puntos de la gráfica (más no que peta)
    logE = zeros(m+1)
    logN = zeros(m+1)

    t1 = dt

    for i in range(0, m+1): 
        aux = len(t1) - 1
        Error, U = ErrorCauchy(F, t1, U0, Esquema, T, N)
        logE[i] = log10(norm(Error[aux, :]))
        logN[i] = log10(float(aux))
        t1 = refine_mesh(t1)

    y = logE[ logE > -12 ] # filtrado de los valores de logE > -12
    x = logN[ 0:len(y) ]
    order = polyfit(x, y, 1)
    print("Order = ", abs(order), "Esquema = ", Esquema)

    return y, x, order

############################################################################
################### REGIONES DE ESTABILIDAD ################################
############################################################################
"""
INPUTS:
    Esquema: Esquema temporal utilizado
    x0: Punto inicial del eje real
    xf: Punto final del eje real
    y0: Punto inicial del eje imaginario
    yf: Punto final del eje imaginario
    N: Número de particiones
OUTPUTS:
    x: equiespaciado de x
    y: equiespaciado de y
    rho: valores de r estables
"""
def Stab(Esquema, x0, xf, y0, yf, N):

    x = linspace(x0, xf, N) # Mallado en el eje real
    y = linspace(y0, yf, N) # Mallado en el eje imaginario

    rho = zeros((N,N)) #

    for i in range(N):
        for j in range(N):
            w = complex(x[i], y[j])
            if Esquema == LeapFrog:
                r = sqrt(Esquema(1, 1, lambda u, t : u*w, 0, 1))
            else:
                r = Esquema(1, 1, lambda u, t : u*w, 0)
            rho[i,j] = abs(r)

    return x, y, rho