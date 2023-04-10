import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gamma, beta
from multiFoxH import compMultiFoxH

from time import time
from sys import exit
from scipy.io import savemat, loadmat

""" implementation of OP equation"""

if __name__ == '__main__':

    L = 2

    """ channel parameters (assuming constant parameters for multiple receivers) """
    alpha = 2   # positive non-linear power parameter
    mu = 2      # number of multipath clusters
    ms = 3      # moderate shadowing
    z_PE = 4 
    psi = mu/(ms-1)
    

    if ms <= 4 / alpha:
        print("condition ms > 4/alpha not met!")
        exit(-1)

    """ SNR threshold """
    gamma_th_dB = 1
    gamma_th = 10 ** (gamma_th_dB / 10)

    """ generate support """
    N = 100
    l_bound_dB = 0
    u_bound_dB = 10
    l_bound = 10 ** (l_bound_dB / 10)
    u_bound = 10 ** (u_bound_dB / 10)
    gamma_bar = np.linspace(l_bound, u_bound, N, dtype=np.float64)
    gamma_bar_dB = 10 * np.log10(gamma_bar)

    """ pre allocation - col vectors """
    OP = np.zeros(N, dtype=np.float64)

    """ PRE COMPUTATIONS - OP analit """
    gamma_coef = (gamma(mu) * gamma(ms)) ** 2

    """ FoxH "arguments" pre computable assuming constant parameters """
    Ai = (1 - ms, 1/alpha)
    Bi = (((z_PE**2/alpha) + 1), 1/alpha)
    Ci = (mu, 1 / alpha)   # a tuple with a pair
    Di = (z_PE**2/alpha, 1/alpha)    # a list of tuples of pairs

    """ compute analytic Outage Probability """

    print("Computing Capacity...")
    OP_analit_start = time()
    for n in range(N):
        start_it = time()
        # a -> top right
        # b -> bot right
        # c -> top left
        # d -> bot left

        """ Parameters for FoxH """
        a = [[(1,1), Ai, Ai, Bi, Bi]] * L
        b = [(Ci, Ci, Di, Di)] * L
        c = [(tuple([0] + [1/2] * L)), (tuple([1] + [-1/2] * L)), (tuple([1] + [-1/2] * L)), (tuple([1] + [1/2] * L))]
        d = [(tuple([1] + [1] * L))]

        mn = [(0, 3)] + [(4, 3)] * L
        pq = [(4, 1)] + [(5, 4)] * L

        Xi = [(psi**(1/alpha))**2 *gamma_bar[n]**(1/2) / (((gamma_bar[n]**2 * (z_PE**2 + 2)** 2))**(1/2))]
        Xis = Xi * L
        z = np.array(np.multiply(Xis, gamma_th))

        param = z, mn, pq, c, d, a, b
        H = np.real(compMultiFoxH(param, nsubdivisions=20, boundaryTol=1e-5))

        OP[n] = 1/(2*np.log(2)) * ( ((z_PE**2)**2) ** L / ((gamma_coef**L) * (alpha**2) ** L )) * H

        end_it = time()
        print("OP[{:>3}] = {:.8e} (L = {:d} => {:.3f} s)".format(n, OP[n], L, end_it - start_it))

    end = time()
    print("OP analit loop time: {:.2f} s" .format(end - OP_analit_start))

    """ save points to .mat file"""
    #save = input("save points? (y/n) ")
    #if save == 'y':

    #filename = 'OP_L{:d}_gamma_bar_{:d}dB.mat'.format(L_MAX, gamma_bar_dB)
    #savemat(filename, dict(L=L_vec, gamma_bar_dB=gamma_bar_dB, gamma_th_dB=gamma_th_dB, OP=OP))
    #print("points saved...")

    """ plotting """
    fig, ax = plt.subplots()
    twin = ax.twinx()

    p_OP, = ax.semilogy(gamma_bar_dB, OP, color='b', linestyle='-')

    ax.set_xlabel("gamma_bar_dB")
    ax.set_ylabel("OP")

    ax.grid()
    #plt.xlim([L_vec[0], L_vec[-1]])
    plt.show()
