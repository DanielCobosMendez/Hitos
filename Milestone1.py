from numpy import array, zeros, linspace, concatenate
from numpy.linalg import norm
import matplotlib.pyplot as plt

# Condiciones iniciales

x0 = 1 # x Inicial
y0 = 0 # y Inicial
dx0 = 0 # Vx inicial
dy0 = 1 # Vy inicial
T = 20 # Periodo
N = 200 # Particiones
dt = T/N # Equiespaciado temporal

# Vector de estado
U0 = array([x0, y0, dx0, dy0]) # Vector de estado inicial

############################################################################
######################## MÉTODO DE EULER IMPLÍCITO #########################
############################################################################

U_Euler = array(zeros((N, len(U0)))) # Vector de estado
U_Euler[0, :] = U0

for i in range(1, N):
    F = concatenate((U_Euler[i-1,2:4], -U_Euler[i-1,0:2]/norm(U_Euler[i-1,0:2])**3)) # Función F
    U_Euler[i,:] = U_Euler[i-1,:] + dt * F # Vector de estado

############################################################################
############################### MÉTODO RK2 #################################
############################################################################

U_RK2 = array(zeros((N, len(U0)))) # Vector de estado
U_RK2[0, :] = U0
K1_RK2 = array(zeros((N, len(U0))))
K2_RK2 = array(zeros((N, len(U0))))

for i in range(1, N):
    K1_RK2[i,:] = concatenate((U_RK2[i-1, 2:4], -U_RK2[i-1, 0:2]/norm(U_RK2[i-1, 0:2])**3))
    K2_RK2[i,:] = concatenate((U_RK2[i-1, 2:4] + dt* K1_RK2[i-1, 2:4], -(U_RK2[i-1, 0:2] + dt* K1_RK2[i-1, 0:2]) / norm(U_RK2[i-1, 0:2] + dt* K1_RK2[i-1, 0:2])**3))
    U_RK2[i,:] = U_RK2[i-1,:] + dt * (K1_RK2[i,:] + K2_RK2[i,:])/ 2

############################################################################
############################### MÉTODO RK4 #################################
############################################################################

U_RK4 = array(zeros((N, len(U0)))) # Vector de estado
U_RK4[0, :] = U0
K1_RK4 = array(zeros((N, len(U0))))
K2_RK4 = array(zeros((N, len(U0))))
K3_RK4 = array(zeros((N, len(U0))))
K4_RK4 = array(zeros((N, len(U0))))

for i in range(1, N):
    K1_RK4[i,:] = concatenate((U_RK4[i-1, 2:4], -U_RK4[i-1, 0:2]/norm(U_RK4[i-1, 0:2])**3))
    K2_RK4[i,:] = concatenate((U_RK4[i-1, 2:4] + dt/2* K1_RK4[i-1, 2:4], -(U_RK4[i-1, 0:2] + dt/2* K1_RK4[i-1, 0:2]) / norm(U_RK4[i-1, 0:2] + dt/2* K1_RK4[i-1, 0:2])**3))
    K3_RK4[i,:] = concatenate((U_RK4[i-1, 2:4] + dt/2* K2_RK4[i-1, 2:4], -(U_RK4[i-1, 0:2] + dt/2* K2_RK4[i-1, 0:2]) / norm(U_RK4[i-1, 0:2] + dt/2* K2_RK4[i-1, 0:2])**3))
    K4_RK4[i,:] = concatenate((U_RK4[i-1, 2:4] + dt* K3_RK4[i-1, 2:4], -(U_RK4[i-1, 0:2] + dt* K3_RK4[i-1, 0:2]) / norm(U_RK4[i-1, 0:2] + dt* K3_RK4[i-1, 0:2])**3))
    U_RK4[i,:] = U_RK4[i-1,:] + dt/6 * (K1_RK4[i,:] + 2 * K2_RK4[i,:] + 2 * K3_RK4[i,:] + K4_RK4[i,:])

############################################################################
############################## GRÁFICAS ####################################
############################################################################

plt.plot(U_Euler[:,0], U_Euler[:,1], label = 'Euler')
plt.plot(U_RK2[:,0], U_RK2[:,1], label = 'RK2')
plt.plot(U_RK4[:,0], U_RK4[:,1], label = 'RK4')
plt.title('Comparación de métodos numéricos')
plt.xlabel('Posición en el eje x')
plt.ylabel('Posición en el eje y')
plt.legend(loc='upper right')
plt.show()

plt.plot(U_Euler[:,2], U_Euler[:,3], label = 'Euler')
plt.plot(U_RK2[:,2], U_RK2[:,3], label = 'RK2')
plt.plot(U_RK4[:,2], U_RK4[:,3], label = 'RK4')
plt.title('Comparación de métodos numéricos')
plt.xlabel('Velocidad en el eje x')
plt.ylabel('Velocidad en el eje y')
plt.legend(loc='upper right')
plt.show()