from numpy import zeros
from ErrorCauchy import ErrorCauchy

def refine_mesh(t1): 
      
      N = len(t1) - 1  
      t2 = zeros( 2*N +1) 

      for i in range(0,N): 
           t2[2*i] =  t1[i]
           t2[2*i+1] = ( t1[i]  + t1[i+1] )/ 2 
      t2[2*N] = t1[N]      

      return t2 

def Conv(F, U0, t, Esquema, T, N):
    """
    Inputs:
        F : Función del problema de Cauchy
        U0 : Vector de estado inicial
        t : Tiempo inicial de partición a integrar
        Esquema: Esquema numérico

    """

    N = len(t)-1
    ptosgraf = 7
    logE = zeros(ptosgraf+1)
    logN = zeros(ptosgraf+1)
    t1 = t

    for i in range(0, ptosgraf+1):
        N = len(t1)-1
        E = ErrorCauchy(F, U0, T, N, Esquema)