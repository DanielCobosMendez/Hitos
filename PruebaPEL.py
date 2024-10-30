from numpy import log, linspace, zeros
import matplotlib.pyplot as plt

def isotermo(x):

    tc = 0.14025*log(1/x)   # x = pc0/pc

    return tc

def adiabatico(x):

    tc = -0.44144*2550532**(-0.4/2.8)*(1-x**(-0.4/2.8))

    return tc



x = linspace(1,1e-15, 500)
tc_iso = zeros(500)
tc_ad = zeros(500)
for i in range(500):

    tc_iso[i] = isotermo(x[i])
    tc_ad[i] = adiabatico(x[i])

plt.axis('equal')
plt.xlabel('pc/pc0')
plt.ylabel('tiempo de cola')
plt.plot(tc_iso, x, '-b', label = 'isotermo')
plt.plot(tc_ad, x, '-r', label = 'adiabático')
plt.show()

# print('el tiempo de cola adiabatico para una x = 0.05 es', adiabatico(0.0000000001))
# print('el tiempo de cola isotermo para una x = 0.05 es', isotermo(0.05))