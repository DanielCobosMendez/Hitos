from numpy import array, zeros

############################################################################
########################### PROBLEMA DE CAUCHY #############################
############################################################################

def Cauchy(F, t, N, U0, Esquema):
    U = array(zeros((N, len(U0)))) # Definición de U
    U[0,:] = U0 # Asignación del vector de estado inicial
    for i in range(1, N):
        U[i, :] = Esquema(U[i-1,:], t/N, F, t) # Llamada al esquema numérico, e integración por cada paso temporal
    return U
