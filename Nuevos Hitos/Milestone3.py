from numpy import concatenate, linspace
from Utiles.Esquemas import ErrorCauchy, Euler, RK2, RK4, EulerI, CrankNicolson, Conv, LeapFrog
from Utiles.Funcion import Funcion, Armonico
import matplotlib.pyplot as plt

############################################################################
##################### CONDICIONES INICIALES ################################
############################################################################
"""
INPUTS:
    FuncionE: Función escogida a resolver 
    T: Periodo
    N: Nº de particiones
OUTPUTS:
    U0: Vector de estado inicial
    dt: malla temporal
"""

FuncionE = Armonico # Función escogida a analizar

if FuncionE == "Armonico":
    U0 = [1, 0]
elif FuncionE == "Funcion":
    r = [1, 0]
    dr = [0, 1]
    U0 = concatenate((r, dr)) # Vector de estado inicial

T = 20 # Periodo
N = 2000 # Particiones
dt = linspace(0, T, N) # Paso temporal

############################################################################
############################# SOLUCIONES ###################################
############################################################################

logE_Euler, logN_Euler, qEuler = Conv(FuncionE, dt, U0, Euler, T, N)
logE_RK2, logN_RK2, qRK2 = Conv(FuncionE, dt, U0, RK2, T, N)
logE_RK4, logN_RK4, qRK4 = Conv(FuncionE, dt, U0, RK4, T, N)
#logE_EulerI, logN_EulerI, qEulerI = Conv(Funcion, dt, U0, EulerI, T, N)
logE_CN, logN_CN, qCN = Conv(FuncionE, dt, U0, CrankNicolson, T, N)
logE_LF, logN_LF, qLF = Conv(FuncionE, dt, U0, LeapFrog, T, N)

############################################################################
############################## GRÁFICAS ####################################
############################################################################

# Gráfica comparación de errores
plt.figure(1)
plt.plot(logN_Euler, logE_Euler, label = 'Euler')
plt.title('Errores')
plt.xlabel('logN')
plt.ylabel('logError')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################

plt.figure(2)
plt.plot(logN_RK2, logE_RK2, label = 'RK2')
plt.title('Errores')
plt.xlabel('logN')
plt.ylabel('logError')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################

plt.figure(3)
plt.plot(logN_RK4, logE_RK4, label = 'RK4')
plt.title('Errores')
plt.xlabel('logN')
plt.ylabel('logError')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################

plt.figure(4)
plt.plot(logN_CN, logE_CN, label = 'CN')
plt.title('Errores')
plt.xlabel('logN')
plt.ylabel('logError')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################

plt.figure(5)
plt.plot(logN_LF, logE_LF, label = 'LF')
plt.title('Errores')
plt.xlabel('logN')
plt.ylabel('logError')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################
plt.show()