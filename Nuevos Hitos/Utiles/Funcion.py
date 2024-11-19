from numpy import concatenate, array
from numpy.linalg import norm

############################################################################
############################ FUNCIÓN ANALIZADA #############################
############################################################################

def Funcion(U, t):
    F = concatenate((U[2:4], -U[0:2]/norm(U[0:2])**3)) # Función F de Euler
    return F


def Armonico(U, t):
    F = array([U[1],-U[0]])
    return F