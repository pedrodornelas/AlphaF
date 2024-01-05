import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gamma, beta
from multiFoxH import compMultiFoxH

from scipy.integrate import simpson

from time import time
from sys import exit

if __name__ == '__main__':

    L = 2

    """ channel parameters (assuming constant parameters for multiple receivers) """

    alpha = 2  # positive non-linear power parameter
    mu = 2  # number of multipath clusters
    ms = 5  # shadowing

    if ms <= 2 / alpha:
        print("condition ms > 2/alpha not met!")
        exit(-1)

    omega = np.sum(np.array([[alpha] * L + [mu] * L]) / 2)

    """ receiver parameters (delta, theta) """
    #delta = 1  # MRC = (1,1); EGC = (1/L, 0.5); SC = (1, +inf)
    #theta = 1000

    """ """
    gamma_th_dB = 10
    gamma_th = 10**(gamma_th_dB/10)

    """ generate support """
    N = 1_00
    l_bound = 1
    u_bound = 100
    gamma_bar = np.linspace(l_bound, u_bound, N)
    gamma_bar_dB = np.multiply(10, np.log10(gamma_bar))

    """ Parameters for FoxH """
    mn = [(0, 1)] + [(1, 2)] * L
    pq = [(1, 1)] + [(2, 2)] * L
    a = [[(1 - ms - mu, 2/alpha), (1 - mu, 2/alpha)]] * L
    b = [[(0, 2/alpha), (-mu, 2/alpha)]] * L
    c = [tuple([-omega] + [1] * L)]
    d = [tuple([1-omega] + [1] * L)]

    """ compute analytic BEP """
    PDF = np.zeros(N, dtype=np.float64)
    start = time()
    for i in range(N):
        Ci = [2/alpha * mu**mu/(gamma(mu) * gamma(ms)) * (1 / ((ms-1)*gamma_bar[i]**(alpha/2)))**mu] * L
        uTheta = [mu ** (2/alpha) / ((ms-1) ** (2/alpha) * gamma_bar[i])] * L
        # print(uTheta)
        z = np.array(np.multiply(uTheta, gamma_th))
        # print(z)
        param = z, mn, pq, c, d, a, b
        H = np.real(compMultiFoxH(param, nsubdivisions=50, boundaryTol=1e-4))
        # print(f'H = {H}')
        PDF[i] = np.prod(np.multiply(Ci, gamma_th**(omega-1) * H))

        # debug
        #print("f[i]     = coef  x prod")
        print("f[{:>4}]  =  {:.2e}".format(i, PDF[i]))
        #print("-----------------------------------------------------------------")

    end = time()
    print("PDF analit loop time: ", end - start)

    area_simpson = simpson(PDF, gamma_bar)
    print("area (simpson): ", area_simpson)

    """ plotting """
    plt.plot(gamma_bar, PDF)
    plt.title("L = " + str(L))
    plt.show()
