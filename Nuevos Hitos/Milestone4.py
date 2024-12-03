from Utiles.Funcion import Armonico
from Utiles.Esquemas import Stab, Euler, RK2, RK4, CrankNicolson, LeapFrog
from numpy import linspace, transpose, shape
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
Funcion = "Armonico" # Función escogida a analizar

estados_iniciales = {"Armonico":[1, 0], "Funcion":[1, 0, 0, 1]}
funciones_iniciales = {"Armonico": Armonico, "Funcion": Funcion}

if Funcion in estados_iniciales:
    U0 = estados_iniciales[Funcion]
    FuncionE = funciones_iniciales[Funcion]
else:
    print("Función no reconocida")
    U0 = None


T = 20 # Periodo
N = 2000 # Particiones
dt = linspace(0, T, N) # Paso temporal

# Puntos iniciales y finales de los ejes real (x) e imaginario (y) de las regiones de estabilidad absoluta
x0 = -4
xf = 4
y0 = -4
yf = 4

############################################################################
############################# SOLUCIONES ###################################
############################################################################

xEuler, yEuler, rhoEuler = Stab(Euler, x0, xf, y0, yf, N)
xRK2, yRK2, rhoRK2 = Stab(RK2, x0, xf, y0, yf, N)
xRK4, yRK4, rhoRK4 = Stab(RK4, x0, xf, y0, yf, N)
xCN, yCN, rhoCN = Stab(CrankNicolson, x0, xf, y0, yf, N)
xLF, yLF, rhoLF = Stab(LeapFrog, x0, xf, y0, yf, N)

############################################################################
############################## GRÁFICAS ####################################
############################################################################

# Crear una figura con 2 filas y 3 columnas
#fig, axs = plt.subplots(2, 3, figsize=(12, 9))

mosaic = """AABBCC
            .DDEE."""

fig, axs = plt.subplot_mosaic(mosaic, figsize=(12, 9), sharey=True)

axs['A'].axis('equal')
# axs[0, 0].axis('equal')
axs['A'].set_title('Región de Estabilidad de Euler')
axs['A'].set_xlabel('Re')
axs['A'].set_ylabel('Im')
axs['A'].contour(xEuler, yEuler, transpose(rhoEuler), linspace(0, 1, 50))
axs['A'].grid()

axs['B'].axis('equal')
#axs[0, 1].axis('equal')
axs['B'].set_title('Región de Estabilidad de RK2')
axs['B'].set_xlabel('Re')
axs['B'].set_ylabel('Im')
axs['B'].contour(xRK2, yRK2, transpose(rhoRK2), linspace(0, 1, 50))
axs['B'].grid()

axs['C'].axis('equal')
#axs[0, 2].axis('equal')
axs['C'].set_title('Región de Estabilidad de RK4')
axs['C'].set_xlabel('Re')
axs['C'].set_ylabel('Im')
axs['C'].contour(xRK4, yRK4, transpose(rhoRK4), linspace(0, 1, 50))
axs['C'].grid()


axs['D'].tick_params('y', labelleft=True)
axs['D'].axis('equal')
#axs[1, 0].axis('equal')
axs['D'].set_title('Región de Estabilidad de CN')
axs['D'].set_xlabel('Re')
axs['D'].set_ylabel('Im')
axs['D'].contour(xCN, yCN, transpose(rhoCN), linspace(0, 1, 50))
axs['D'].grid()

axs['E'].axis('equal')
#axs[1, 1].axis('equal')
axs['E'].set_title('Región de Estabilidad de Leap Frog')
axs['E'].set_xlabel('Re')
axs['E'].set_ylabel('Im')
axs['E'].contour(xLF, yLF, transpose(rhoLF), linspace(0, 1, 50))
axs['E'].grid()

plt.tight_layout()

#axs[1, 2].axis('off') # Oculta el subplot vacío

plt.show()