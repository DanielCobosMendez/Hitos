from numpy import array, zeros, linspace
from TemporalScheme.Esquemas import Euler, RK2, RK4, EulerI, CrankNicholson
from TemporalScheme.Cauchy import Cauchy
from TemporalScheme.Funcion import Funcion
from TemporalScheme.Error import Error
from TemporalScheme.ErrorCauchy import ErrorCauchy
from TemporalScheme.Convergencia import Conv
import matplotlib.pyplot as plt

############################################################################
######################## CONDICIONES INICIALES #############################
############################################################################
#         ->->->->-> Lo único que hay que cambiar <-<-<-<-

T = 0.5 # Periodo
N = 10 # Nº de particiones de U1



x0 = 1 # x Inicial
y0 = 0 # y Inicial
dx0 = 0 # Vx inicial
dy0 = 1 # Vy inicial
U0 = array([x0, y0, dx0, dy0]) # Vector de estado inicial euler
#U0 = array([x0, dx0])


############################################################################
############################# SOLUCIONES ###################################
############################################################################

E_Euler = ErrorCauchy(Funcion, U0, T, N, Euler)
E_RK2 = ErrorCauchy(Funcion, U0, T, N, RK2)
E_RK4 = ErrorCauchy(Funcion, U0, T, N, RK4)
E_CN = ErrorCauchy(Funcion, U0, T, N, CrankNicholson)

print(E_Euler)
############################################################################
############################ CONVERGENCIA ##################################
############################################################################

q_Euler, logNEuler = Conv(Funcion, U0, T, N, Euler)
#q_RK2 = Conv(Funcion, U0, T, N, RK2)
#q_RK4 = Conv(Funcion, U0, T, N, RK4)
#q_CN = Conv(Funcion, U0, T, N, CrankNicholson)

############################################################################
############################## GRÁFICAS ####################################
############################################################################

# Gráfica comparación de posiciones
plt.plot(linspace(0, T, N), E_RK2[:,0], label = 'RK2')
plt.plot(linspace(0, T, N), E_RK4[:,0], label = 'RK4')
plt.plot(linspace(0, T, N), E_CN[:, 0], label = 'Crank-Nicholson')
plt.title('Comparación error de métodos numéricos')
plt.xlabel('N')
plt.ylabel('Error')
plt.legend(loc='upper right')
plt.grid() #############
#plt.axis('equal') ###################
plt.ylim(-0.003,0.003)
plt.show()

# Gráfica comparación de posiciones
plt.plot(linspace(0, T, N), E_Euler[:,0], label = 'Euler')
#plt.plot(E_EulerI[:,0], N, label = 'Euler Inverso')
plt.title('Comparación error de métodos numéricos')
plt.xlabel('N')
plt.ylabel('Error')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################
plt.show()

