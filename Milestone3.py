from numpy import array, zeros, linspace
from TemporalScheme.Esquemas import Euler, RK2, RK4, EulerI, CrankNicholson
from TemporalScheme.Cauchy import Cauchy
from TemporalScheme.Funcion import Funcion
from TemporalScheme.Error import Error
import matplotlib.pyplot as plt

############################################################################
######################## CONDICIONES INICIALES #############################
############################################################################
#         ->->->->-> Lo único que hay que cambiar <-<-<-<-

T = 20 # Periodo
N2 = 2000 # Nº de particiones de U1
N = N2*2 # Nº de particiones de U2


x0 = 1 # x Inicial
y0 = 0 # y Inicial
dx0 = 0 # Vx inicial
dy0 = 1 # Vy inicial
U0 = array([x0, y0, dx0, dy0]) # Vector de estado inicial

############################################################################
############################# SOLUCIONES ###################################
############################################################################

U_Euler_1 = Cauchy(Funcion, T, N2, U0, Euler)
U_RK2_1 = Cauchy(Funcion, T, N2, U0, RK2)
U_RK4_1 = Cauchy(Funcion, T, N2, U0, RK4)
#U_EulerI_1 = Cauchy(Funcion, T, N2, U0, EulerI)
U_CN_1 = Cauchy(Funcion, T, N2, U0, CrankNicholson)

U_Euler_2 = Cauchy(Funcion, T, N, U0, Euler)
U_RK2_2 = Cauchy(Funcion, T, N, U0, RK2)
U_RK4_2 = Cauchy(Funcion, T, N, U0, RK4)
#U_EulerI_2 = Cauchy(Funcion, T, N, U0, EulerI)
U_CN_2 = Cauchy(Funcion, T, N, U0, CrankNicholson)



############################################################################
############################# SOLUCIONES ###################################
############################################################################

E_Euler = Error(U_Euler_1, U_Euler_2, U0, N2, Euler)
E_RK2 = Error(U_RK2_1, U_RK2_2, U0, N2, RK2)
E_RK4 = Error(U_RK4_1, U_RK4_2, U0, N2, RK4)
#U_EulerI = Error(U_Euler_1, U_Euler_2, U0, N2, Euler)
E_CN = Error(U_CN_1, U_CN_2, U0, N2, CrankNicholson)

############################################################################
############################## GRÁFICAS ####################################
############################################################################

# Gráfica comparación de posiciones
#plt.plot(linspace(0, T, N2), E_Euler[:,0], label = 'Euler')
#plt.plot(E_EulerI[:,0], N, label = 'Euler Inverso')
plt.plot(linspace(0, T, N2), E_RK2[:,0], label = 'RK2')
plt.plot(linspace(0, T, N2), E_RK4[:,0], label = 'RK4')
plt.plot(linspace(0, T, N2), E_CN[:, 0], label = 'Crank-Nicholson')
plt.title('Comparación error de métodos numéricos')
plt.xlabel('N')
plt.ylabel('Error')
plt.legend(loc='upper right')
plt.grid() #############
#plt.axis('equal') ###################
plt.ylim(-0.003,0.003)
plt.show()

# Gráfica comparación de posiciones
plt.plot(linspace(0, T, N2), E_Euler[:,0], label = 'Euler')
#plt.plot(E_EulerI[:,0], N, label = 'Euler Inverso')
plt.title('Comparación error de métodos numéricos')
plt.xlabel('N')
plt.ylabel('Error')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################
plt.show()