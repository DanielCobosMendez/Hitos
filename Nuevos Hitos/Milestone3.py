from numpy import concatenate, linspace
from Utiles.Esquemas import ErrorCauchy, Euler, RK2, RK4, EulerI, CrankNicolson, Conv
from Utiles.Funcion import Funcion
import matplotlib.pyplot as plt

# Condiciones iniciales Euler

r = [1, 0]
dr = [0, 1]
U0 = concatenate((r, dr)) # Vector de estado inicial

T = 20 # Periodo
N = 2000 # Particiones
dt = linspace(0, T, N) # Paso temporal

############################################################################
############################# SOLUCIONES ###################################
############################################################################

logE_Euler, logN_Euler, qEuler = Conv(Funcion, dt, U0, Euler, T, N)
logE_RK2, logN_RK2, qRK2 = Conv(Funcion, dt, U0, RK2, T, N)
logE_RK4, logN_RK4, qRK4 = Conv(Funcion, dt, U0, RK4, T, N)
#logE_EulerI, logN_EulerI, qEulerI = Conv(Funcion, dt, U0, EulerI, T, N)
logE_CN, logN_CN, qCN = Conv(Funcion, dt, U0, CrankNicolson, T, N)

############################################################################
############################## GRÁFICAS ####################################
############################################################################

# Gráfica comparación de errores
plt.plot(logN_Euler, logE_Euler, label = 'Euler')
plt.plot(logN_RK2, logE_RK2, label = 'RK2')
plt.plot(logN_RK4, logE_RK4, label = 'RK4')
plt.plot(logN_CN, logE_CN, label = 'CN')
plt.title('Comparación de erroes')
plt.xlabel('logN')
plt.ylabel('logError')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################
plt.show()