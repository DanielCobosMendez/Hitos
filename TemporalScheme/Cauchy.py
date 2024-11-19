from numpy import array, zeros

############################################################################
########################### PROBLEMA DE CAUCHY #############################
############################################################################

def Cauchy(F, t, U0, Esquema, q = None, Tol = None ):

    """
    INPUTS:
        F(U, t): Función dU/dt = F(U, t)
        t: Partición temporal t (vector de longitud N+1)
        U0: Vector de estado inicial
        Esquema: Esquema temporal
        q: Orden (opcional)
        Tol: Tolerancia (opcional)
    OUTPUTS:
        U: matriz[N+1, Nv] Son los Nv valores en el paso temporal N
    """

    N, Nv =  len(t)-1, U0.size
    U = zeros( (N+1, Nv), dtype=type(U0) )

    U[0,:] = U0

    for n in range(N):

        if q != None and Tol !=None: 
          U[n+1,:] = Esquema( U[n, :], t[n+1] - t[n], t[n],  F, q, Tol ) 

        else: 
           U[n+1,:] = Esquema( U[n, :], t[n+1] - t[n], t[n],  F ) 

    return U

    #U = array(zeros((N, len(U0)))) # Definición de U
    #U[0,:] = U0 # Asignación del vector de estado inicial
    #for i in range(1, N):
    #    U[i, :] = Esquema(U[i-1,:], t/N, F, t) # Llamada al esquema numérico, e integración por cada paso temporal
    #return U
#