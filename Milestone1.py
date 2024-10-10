from numpy import array, zeros, linspace
import matplotlib.pyplot as plt

# Método de Euler implícito

U = array([1, 0, 0, 1]) # Vector de estado de condiciones iniciales

N = 200 # Nº de particiones
T = 20 # Periodo

x = array(zeros(N)) # Array de los valores de x
y = array(zeros(N)) # Array de los valores de y
dx = array(zeros(N)) # Array de los valores de x punto
dy = array(zeros(N)) # Array de los valores de y punto
x[0] = U[0]
y[0] = U[1]
dx[0] = U[2]
dy[0] = U[3]

for i in range(1, N):
    F = array([U[2], U[3], -U[0]/(U[0]**2+U[1]**2)**1.5, -U[1]/(U[0]**2+U[1]**2)**1.5]) # Vector F, derivada del vector de estado
    U = U + T/N * F
    x[i] = U[0]
    y[i] = U[1]
    dx[i] = U[2]
    dy[i] = U[3]

plt.plot(x, y)
plt.title('Órbita de Kepler según el método de Euler')
plt.xlabel('Eje x')
plt.ylabel('Eje y')
plt.show()

plt.plot(dx, dy)
plt.title('Orbita de Kepler según el método de Euler')
plt.xlabel('Velocidad según el eje x')
plt.ylabel('Velocidad según el eje y')
plt.show()