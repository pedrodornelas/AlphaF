import matplotlib.pyplot as plt
import numpy as np
from scipy.special import gamma, beta
from multiFoxH import compMultiFoxH

import os
from time import time
from sys import exit
from scipy.io import savemat, loadmat

def product(b, Ui, j):
    #print(b[0])
    #print(b[0][j])
    prod = 1
    for i in range(len(b[0])):
        if i != j:
            prod *= gamma(b[0][i][0]-b[0][i][1]*Ui)
            # print('gamma('+str(b[0][i][0])+'-'+str(b[0][i][1])+'*'+str(Ui)+')')
    return prod

""" implementation of OP equation"""

if __name__ == '__main__':

    """ plotting """
    fig, ax = plt.subplots()
    # twin = ax.twinx()

    L = [1,2,3,4]
    # L = [1,2]

    for k in range(len(L)):
        """ channel parameters (assuming constant parameters for multiple receivers) """
        alpha = [1.5,2.3]   # positive non-linear power parameter
        mu = [2,2]          # number of multipath clusters
        ms = [3,4]          # moderate shadowing
        z_PE = [0.8,1.5]    # pointing error
        """ SNR threshold """
        gamma_th_dB = 5
        gamma_th = 10 ** (gamma_th_dB / 10)


        # psi = mu/(ms-1)
        psi1 = mu[0]/(ms[0]-1)
        psi2 = mu[1]/(ms[1]-1)

        for i in range(len(ms)):
            if ms[i] <= 4/alpha[i]:
                print("condition ms > 4/alpha not met!")
                exit(-1)

        """ generate support """
        N = 100
        l_bound_dB = 0
        u_bound_dB = 30
        l_bound = 10 ** (l_bound_dB / 10)
        u_bound = 10 ** (u_bound_dB / 10)
        # gamma_bar = np.linspace(l_bound, u_bound, N, dtype=np.float64)
        gamma_bar_dB = np.linspace(l_bound_dB, u_bound_dB, N, dtype=np.float64)
        #gamma_bar_dB = 10 * np.log10(gamma_bar)
        gamma_bar = 10 ** (gamma_bar_dB / 10)

        """ pre allocation - col vectors """
        OP = np.zeros(N, dtype=np.float64)
        OP_asymptotic = np.zeros(N, dtype=np.float64)

        """ FoxH "arguments" pre computable assuming constant parameters """
        # Ai = (1 - ms, 1/alpha)
        Ai1 = (1-ms[0], 1/alpha[0])
        Ai2 = (1-ms[1], 1/alpha[1])
        
        # Bi = (((z_PE**2/alpha) + 1), 1/alpha)
        Bi1 = (((z_PE[0]**2/alpha[0]) + 1), 1/alpha[0])
        Bi2 = (((z_PE[1]**2/alpha[1]) + 1), 1/alpha[1])
        
        # Ci = (mu, 1 / alpha)   # a tuple with a pair 
        Ci1 = (mu[0], 1 / alpha[0])   # a tuple with a pair 
        Ci2 = (mu[1], 1 / alpha[1])   # a tuple with a pair 
        
        # Di = (z_PE**2/alpha, 1/alpha)    # a list of tuples of pairs
        Di1 = (z_PE[0]**2/alpha[0], 1/alpha[0])    # a list of tuples of pairs
        Di2 = (z_PE[1]**2/alpha[1], 1/alpha[1])    # a list of tuples of pairs

        Epsilon1 = (1, 1)
        Epsilon2 = (0, 1/2)
        Epsilon3 = (1, 1/2)

        """ compute analytic Outage Probability """
        # a -> top right
        # b -> bot right
        # c -> top left
        # d -> bot left

        """ Parameters for FoxH """
        a = [[(1,1), Ai1, Ai2, Bi1, Bi2]] * L[k]
        b = [(Ci1, Ci2, Di1, Di2)] * L[k]
        c = [tuple([1] + [1/2] * L[k])]                                 # Epsilon3
        d = [(tuple([1] + [1] * L[k])), (tuple([0] + [1/2] * L[k]))]    # Epsilon1 and Epsilon2

        mn = [(0, 1)] + [(4, 3)] * L[k]
        pq = [(1, 2)] + [(5, 4)] * L[k]

        arr_Ui = []
        for tupla in b[0]:
            arr_Ui.append(tupla[0]/tupla[1])

        Ui = min(arr_Ui)
        j = arr_Ui.index(Ui) #minimum position

        
        # print('Arr_UI: '+str(arr_Ui))
        # print('Min: '+str(Ui))
        # print('Pos min: '+str(j))

        # b[0][j][0]
        # print(gamma(b[0][j][0]-b[0][j][1]*Ui))

        Beta = b[0][j][1]
        # print('Beta: '+str(Beta))

        """ PRE COMPUTATIONS - OP analit """
        gamma_coef1 = gamma(mu[0]) * gamma(ms[0])
        gamma_coef2 = gamma(mu[1]) * gamma(ms[1])
        preH = 1/2 * (z_PE[0]*z_PE[1])**L[k] / (alpha[0]*alpha[1]*gamma_coef1*gamma_coef2)**L[k]

        print("Computing OP...")
        OP_analit_start = time()
        for n in range(N):
            start_op = time()

            Xi = [(psi1**(1/alpha[0]))*(psi2**(1/alpha[1])) / ((gamma_bar[n]**2 * (z_PE[0]**2 + 2) * (z_PE[1]**2 + 2))**(1/2))]
            Xis = Xi * L[k]
            # print(Xis)
            z = np.array(np.multiply(Xis, gamma_th**(1/2)))

            param = z, mn, pq, c, d, a, b
            H = np.real(compMultiFoxH(param, nsubdivisions=35, boundaryTol=1e-5))

            OP[n] = preH * H
            end_op = time()

            start_asymp = time()

            Psi_Ui = gamma(1-Epsilon3[0]+Epsilon3[1]*L[k]*Ui)/(gamma(1-Epsilon2[0]+Epsilon2[1]*L[k]*Ui)*gamma(1-Epsilon1[0]+Epsilon1[1]*L[k]*Ui))
            Phi_Ui = product(b, Ui, j) * gamma(1-Ai1[0]+Ai1[1]*Ui) * gamma(1-Ai2[0]+Ai2[1]*Ui) * gamma(Ui) / (gamma(Bi1[0]-Bi1[1]*Ui)*gamma(Bi2[0]-Bi2[1]*Ui))
            OP_asymptotic[n] = preH * Psi_Ui * (Phi_Ui**L[k]) * ((Xi[0]*(gamma_th**(1/2)))**(Ui*L[k])) / (Beta**L[k])

            end_asymp = time()

            print("OP [{:>3}] = {:.8e} (L = {:d} => {:.3f} s)".format(n, OP[n], L[k], end_op - start_op))
            print("ASY[{:>3}] = {:.8e} (L = {:d} => {:.3f} s)".format(n, OP_asymptotic[n], L[k], end_asymp - start_asymp))

        end = time()
        print("OP analit loop time: {:.2f} s" .format(end - OP_analit_start))

        """ save points to .mat file"""
        #save = input("save points? (y/n) ")
        #if save == 'y':

        filename = 'MAT/OP_L_{:d}_gamma_th_{:d}_dB.mat'.format(L[k], gamma_th_dB)
        savemat(filename, dict(L=L, gamma_bar_dB=gamma_bar_dB, gamma_th_dB=gamma_th_dB, OP=OP, OP_Asymp=OP_asymptotic))
        print("points saved...")

        color = ['b','g','r','m']

        p_OP = ax.semilogy(gamma_bar_dB, OP, color=color[k], linestyle='-', label=('L = '+str(L[k])))
        if k == len(L)-1:
            p_OP_asymp = ax.semilogy(gamma_bar_dB, OP_asymptotic, color='k', linestyle='--', label='Asymptotic')
        else:
            p_OP_asymp = ax.semilogy(gamma_bar_dB, OP_asymptotic, color='k', linestyle='--')

    lengend = ax.legend(loc='lower left')

    ax.set_xlabel("gamma_bar_dB")
    ax.set_ylabel("OP")

    ax.grid()
    
    # plt.ylim([10**(-5), 10**0])
    plt.show()
