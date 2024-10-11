from numpy import array, zeros, linspace, concatenate
from numpy.linalg import norm
import matplotlib.pyplot as plt

############################################################################
############################ FUNCIÓN ANALIZADA #############################
############################################################################

def Funcion(U, t):
    F = concatenate((U[2:4], -U[0:2]/norm(U[0:2])**3)) # Función F de Euler
    return F

############################################################################
######################## MÉTODO DE EULER IMPLÍCITO #########################
############################################################################

def Euler(U, dt, F, t): # U(n+1) = U(n) + dt * F

    return U + dt * F(U, t)

############################################################################
############################### MÉTODO RK2 #################################
############################################################################

def RK2(U, dt, F, t): # U(n+1) = U(n) + dt/2 * (K1 + K2)
    K1_RK2 = F(U, t)
    K2_RK2 = F(U + K1_RK2*dt, t + dt)
    return U + dt/2 * (K1_RK2 + K2_RK2)

############################################################################
############################### MÉTODO RK4 #################################
############################################################################

def RK4(U, dt, F, t): # U(n+1) = U(n) + dt/6 * (K1 + 2*K2 + 2*K3 + K4)
    K1_RK4 = F(U, t)
    K2_RK4 = F(U + K1_RK4*dt/2, t + dt/2)
    K3_RK4 = F(U + K2_RK4*dt/2, t + dt/2)
    K4_RK4 = F(U + K3_RK4*dt, t + dt)
    return U + dt/6 * (K1_RK4 + 2 * K2_RK4 + 2 * K3_RK4 + K4_RK4)

############################################################################
########################### PROBLEMA DE CAUCHY #############################
############################################################################

def Cauchy(F, t, N, U0, Esquema):
    U = array(zeros((N, len(U0)))) # Definición de U
    U[0,:] = U0 # Asignación del vector de estado inicial
    for i in range(1, N):
        U[i, :] = Esquema(U[i-1,:], t/N, F, t) # Llamada al esquema numérico, e integración por cada paso temporal
    return U

############################################################################
######################## CONDICIONES INICIALES #############################
############################################################################
#         ->->->->-> Lo único que hay que cambiar <-<-<-<-

T = 20 # Periodo
N = 200 # Nº de particiones

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

############################################################################
############################## GRÁFICAS ####################################
############################################################################

# Gráfica comparación de posiciones
plt.plot(U_Euler[:,0], U_Euler[:,1], label = 'Euler')
plt.plot(U_RK2[:,0], U_RK2[:,1], label = 'RK2')
plt.plot(U_RK4[:,0], U_RK4[:,1], label = 'RK4')
plt.title('Comparación de métodos numéricos')
plt.xlabel('Posición en el eje x')
plt.ylabel('Posición en el eje y')
plt.legend(loc='upper right')
plt.show()

# Gráfica comparación perfil de velocidades
plt.plot(U_Euler[:,2], U_Euler[:,3], label = 'Euler')
plt.plot(U_RK2[:,2], U_RK2[:,3], label = 'RK2')
plt.plot(U_RK4[:,2], U_RK4[:,3], label = 'RK4')
plt.title('Comparación de métodos numéricos')
plt.xlabel('Velocidad en el eje x')
plt.ylabel('Velocidad en el eje y')
plt.legend(loc='upper right')
plt.show()

plt.plot(linspace(0, T, N), U_Euler[:,2], label = 'Vx Euler')
plt.plot(linspace(0, T, N), U_RK2[:,2], label = 'Vx RK2')
plt.plot(linspace(0, T, N), U_RK4[:,2], label = 'Vx RK4')
plt.title('Comparación de velocidades en X')
plt.xlabel('Tiempo')
plt.ylabel('Velocidad')
plt.legend(loc='upper right')
plt.show()

plt.plot(linspace(0, T, N), U_Euler[:,3], label = 'Vy Euler')
plt.plot(linspace(0, T, N), U_RK2[:,3], label = 'Vy RK2')
plt.plot(linspace(0, T, N), U_RK4[:,3], label = 'Vy RK4')
plt.title('Comparación de velocidades en y')
plt.xlabel('Tiempo')
plt.ylabel('Velocidad')
plt.legend(loc='upper right')
plt.show()