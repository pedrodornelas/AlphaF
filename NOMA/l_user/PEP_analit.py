import numpy as np
from scipy.special import gamma, binom
from multiFoxH import parseArgsToMultiH
import math

pi = math.pi

def PEP_analit(L: int, l: int, params: list, gammaBar: np.array) -> list[float]:
    points = len(gammaBar)
    # print(points)

    alpha = params[0]
    mu = params[1]
    ms = params[2]
    z = params[3]
    theta = 1

    gamma_ratio = (gamma(L+1)) / (gamma(l)*gamma(L-l+1))
    preH = gamma_ratio * (alpha / (4*(pi**(1/2))))
    
    sum = 0

    for k in range(L-l):
        # print(k)
        binomial = binom((L-l), k)
        signal = (-1)**k
        frac = ((z**2) / (alpha*gamma(mu)*gamma(ms)))**(l+k)

        Psi = (mu/(ms-1)) ** (1/alpha)
        Xi = (Psi*( (z*(2**(1/2)))**alpha )) / (((theta**2) * gammaBar * (z**2 + 2))**(alpha/2))
        # print(Xi)

        H = parseArgsToMultiH(params, l, k, Xi, points)
        # print(len(H))

        sum = sum + (binomial*signal*frac*H)

    PEP = preH * sum
    # print(len(PEP))

    return PEP