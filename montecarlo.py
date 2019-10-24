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
        print(c)

    sigma = np.std(yCoords) 
    print(sigma)

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
    plotCoords(0, 2, 20, 10000000, [mc.nfEquation])





