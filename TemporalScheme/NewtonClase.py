from numpy import array, exp, linspace

def Newton(F, x_0, Fprima = None, tol = 1e-4, maxiter=50):
    def Fp(x):
        if Fprima==None:
            delta = 1e-4
            return (F(x+delta)-F(x-delta))/(2*delta)
        else:
            return Fprima(x)
    
    xn = x_0
    Error = tol + 1
    iter = 0

    while Error > tol and iter<maxiter:

        xn1 = xn - F(xn)/Fp(xn)
        Error = abs(xn-xn1)
        xn = xn1

        iter += 1 # Sumatorio
        print('Error:', Error)
    print('Número de iteraciones:', iter)
    return xn

def jorge(x):
    return exp(x)-2*x-2 # Al usar la exp de numpy, se aceptan inputs tanto como escalar o vectorial (lista)


def dif_jorge(x):
    return exp(x)-2

x = linspace(-2, 2, 1000)
y = jorge(x)

Sol = Newton(jorge, -10, dif_jorge)
print(Sol)
print("Residual=",jorge(Sol))

Sol2 = Newton(jorge, -10)
print(Sol2)
print("Residual=",jorge(Sol2))