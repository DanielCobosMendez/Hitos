from numpy import array, zeros, linspace, concatenate
from numpy.linalg import norm
import matplotlib.pyplot as plt
from scipy.optimize import newton
from TemporalScheme.Esquemas import Euler, RK2, RK4, EulerI, CrankNicholson
from TemporalScheme.Cauchy import Cauchy
from TemporalScheme.Funcion import Funcion

############################################################################
######################## CONDICIONES INICIALES #############################
############################################################################
#         ->->->->-> Lo único que hay que cambiar <-<-<-<-

T = 20 # Periodo
N = 20000 # Nº de particiones

x0 = 1 # x Inicial
y0 = 0 # y Inicial
dx0 = 0 # Vx inicial
dy0 = 1 # Vy inicial
U0 = array([x0, y0, dx0, dy0]) # Vector de estado inicial

############################################################################
############################# SOLUCIONES ###################################
############################################################################

U_Euler = Cauchy(Funcion, T, N, U0, Euler)
U_RK2 = Cauchy(Funcion, T, N, U0, RK2)
U_RK4 = Cauchy(Funcion, T, N, U0, RK4)
U_EulerI = Cauchy(Funcion, T, N, U0, EulerI)
U_CN = Cauchy(Funcion, T, N, U0, CrankNicholson)

############################################################################
############################## GRÁFICAS ####################################
############################################################################

# Gráfica comparación de posiciones
plt.plot(U_Euler[:,0], U_Euler[:,1], label = 'Euler')
plt.plot(U_EulerI[:,0], U_EulerI[:,1], label = 'Euler Inverso')
plt.plot(U_RK2[:,0], U_RK2[:,1], label = 'RK2')
plt.plot(U_RK4[:,0], U_RK4[:,1], label = 'RK4')
plt.plot(U_CN[:, 0], U_CN[:, 1], label = 'Crank-Nicholson')
plt.title('Comparación de métodos numéricos')
plt.xlabel('Posición en el eje x')
plt.ylabel('Posición en el eje y')
plt.legend(loc='upper right')
plt.grid() #############
plt.axis('equal') ###################
plt.show()

# Gráfica comparación perfil de velocidades
plt.plot(U_Euler[:,2], U_Euler[:,3], label = 'Euler')
plt.plot(U_EulerI[:,2], U_EulerI[:,3], label = 'Euler Inverso')
plt.plot(U_RK2[:,2], U_RK2[:,3], label = 'RK2')
plt.plot(U_RK4[:,2], U_RK4[:,3], label = 'RK4')
plt.plot(U_CN[:, 2], U_CN[:, 3], label = 'Crank-Nicholson')
plt.title('Comparación de métodos numéricos')
plt.xlabel('Velocidad en el eje x')
plt.ylabel('Velocidad en el eje y')
plt.legend(loc='upper right')
plt.show()

plt.plot(linspace(0, T, N), U_Euler[:,0], label = 'x Euler')
plt.plot(linspace(0, T, N), U_EulerI[:,0], label = 'x Euler Inverso')
plt.plot(linspace(0, T, N), U_RK2[:,0], label = 'x RK2')
plt.plot(linspace(0, T, N), U_RK4[:,0], label = 'x RK4')
plt.plot(linspace(0, T, N), U_CN[:, 0], label = 'Crank-Nicholson')
plt.title('Comparación de posiciones en X')
plt.xlabel('Tiempo')
plt.ylabel('Velocidad')
plt.legend(loc='upper right')
plt.show()

plt.plot(linspace(0, T, N), U_Euler[:,1], label = 'y Euler')
plt.plot(linspace(0, T, N), U_EulerI[:,1], label = 'y Euler Inverso')
plt.plot(linspace(0, T, N), U_RK2[:,1], label = 'y RK2')
plt.plot(linspace(0, T, N), U_RK4[:,1], label = 'y RK4')
plt.plot(linspace(0, T, N), U_CN[:, 1], label = 'Crank-Nicholson')
plt.title('Comparación de posiciones en y')
plt.xlabel('Tiempo')
plt.ylabel('Velocidad')
plt.legend(loc='upper right')
plt.show()

plt.plot(linspace(0, T, N), U_Euler[:,2], label = 'Vx Euler')
plt.plot(linspace(0, T, N), U_EulerI[:,2], label = 'Vx Euler Inverso')
plt.plot(linspace(0, T, N), U_RK2[:,2], label = 'Vx RK2')
plt.plot(linspace(0, T, N), U_RK4[:,2], label = 'Vx RK4')
plt.plot(linspace(0, T, N), U_CN[:, 2], label = 'Crank-Nicholson')
plt.title('Comparación de velocidades en X')
plt.xlabel('Tiempo')
plt.ylabel('Velocidad')
plt.legend(loc='upper right')
plt.show()

plt.plot(linspace(0, T, N), U_Euler[:,3], label = 'Vy Euler')
plt.plot(linspace(0, T, N), U_EulerI[:,3], label = 'Vy Euler Inverso')
plt.plot(linspace(0, T, N), U_RK2[:,3], label = 'Vy RK2')
plt.plot(linspace(0, T, N), U_RK4[:,3], label = 'Vy RK4')
plt.plot(linspace(0, T, N), U_CN[:, 3], label = 'Crank-Nicholson')
plt.title('Comparación de velocidades en y')
plt.xlabel('Tiempo')
plt.ylabel('Velocidad')
plt.legend(loc='upper right')
plt.show()