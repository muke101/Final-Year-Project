import random
import numpy as np
from matplotlib import pyplot as plt

def CaEquation(x):
    return x**2

def nfEquation(x, k, z, phi):
    fc = z*(1+2*((1-z)/z)**(1/2)*np.tan(k)*np.cos(phi)+np.tan(k)**2*(1-z)/z)**(1-x/2)+(1-z)*(1-2*(z/(1-z))**(1/2)*np.tan(k)*np.cos(phi)+(z*np.tan(k)**2)/(1-z))**(1-x/2) 
    Hq = 1-((z*(1-z))/(1+np.tan(k)**2))*(2*np.cos(phi)+(np.tan(k)*(1-2*z))/((z*(1-z))**(1/2)))**2
    return (np.log(fc)*Hq)/(4*np.pi)
    return (np.log(fc)*Hq)/(4*np.pi)

def monteCarlo(equation, x, N):
    I = 0
    V = np.pi**2
    c = 0
    results = np.empty(N-1)
    for i in range(1,N):
        phi = random.uniform(0, 2*np.pi)
        z = random.uniform(0,1)
        k = random.uniform(0,np.pi/2)
        results[c] = equation(x, k, z, phi)
        print(I)
        I += results[c]
        c += 1
    
    integral = (V*I)/N
    u = np.mean(results)
    
    variance = 0
    for i in results:
        variance += (i-u)**2
    sigma = (variance/N)**(1/2)
    
    return integral

def yCoords(equation, xCoords):
    yCoords = []
    for x in xCoords:
        yCoords.append(monteCarlo(equation, 0.8, 10000))

    return yCoords

def plotCoords(x1, x2, n, equations):
    xCoords = np.linspace(x1, x2, n)

    for eq in equations:
        plt.plot(xCoords, yCoords(eq, xCoords))

    #plt.show()

#plotCoords(0, 2, 20, [nfEquation])
print(nfEquation(0.8, np.pi/4, 0.5, np.pi))





