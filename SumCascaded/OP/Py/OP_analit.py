import numpy as np
from scipy.special import gamma
from multiFoxH import parseArgsFromMatlab

def OP_analit(L: int, N: int, params: list, gamma_th: float, gammaBar: np.array):
    points = len(gammaBar)
    # print(points)

    Xi = 1
    preH = 1

    for c in range(N):
        alpha = params[c][0]
        mu = params[c][1]
        ms = params[c][2]
        z = params[c][3]

        Psi = (mu/(ms-1)) ** (1/alpha)
        Xi = Xi * ((Psi * z) / ((gammaBar[:, c] * (z**2+2))**(1/2)))

        preH = preH * (z**2 / (alpha * gamma(mu) * gamma(ms)))

    Xi = Xi * (gamma_th**(1/2))
    preH = (preH ** L) / 2
    # print(preH)

    H = parseArgsFromMatlab(params, N, L, Xi)
    # print(len(H))
    # print(H)

    OP = preH * H
    # print(len(OP))
    # print(OP)

    return OP