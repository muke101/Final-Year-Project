import numpy as np 
import random

def equation(phi_k, z, k, x):
    while k == 1:
        k = random.uniform(0,1)
    u2 = k/(1-k)
    phi = phi_k*np.pi
    u = u2**0.5

    fc = z*(1+2*((1-z)/z)**0.5*u*np.cos(phi)+(u2*(1-z))/z)**(1-x/2)+(1-z)*(1-2*(z/(1-z))**0.5*u*np.cos(phi)+z/(1-z)*u2)**(1-x/2)
    Hq = 1-(z*(1-z)/(1+u2))*(2*np.cos(phi)+((1-2*z)*u)/(z*(1-z))**0.5)**2

    return 1/(u2*(1+u2))*1/2*1/(2*np.pi)*Hq*np.log(fc)

def monteCarlo(equation, x, N):
    I = 0
    
    for i in range(N):
        phi_k = random.uniform(0,1)
        z = random.uniform(0,1)
        k = random.uniform(0,1)
        I += equation(phi_k, z, k, x)

    return I/N

print(monteCarlo(equation, 0, 10000))

