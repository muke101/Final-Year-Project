import random
import ctypes
import numpy as np
from matplotlib import pyplot as plt

def yCoords(equation, xCoords, N):
    yCoords = np.empty(len(xCoords))

    N = ctypes.c_ulonglong(N)

    for c, x in enumerate(xCoords):
        x = ctypes.c_double(x)
        yCoords[c] = mc.monteCarlo(equation, x, N)
        output.write(str(c+1)+": "+str(yCoords[c]))
        output.write('\n')
        print(c+1)

    sigma = np.std(yCoords) 

    output.write("std: "+str(sigma))
    output.write('\n')

    return yCoords

def plotCoords(x1, x2, n, N, equations):
    xCoords = np.linspace(x1, x2, n)

    for eq in equations:
        plt.plot(xCoords, yCoords(eq, xCoords, N))

    plt.show()

if __name__ == '__main__':
    ctypes.cdll.LoadLibrary("./montecarlo.so")
    mc = ctypes.CDLL("./montecarlo.so")
    mc.monteCarlo.restype = ctypes.c_double
    output = open('results', 'w')
    plotCoords(0, 2, 50, 10000000, [mc.nfEquation])
    output.close()





