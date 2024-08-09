import numpy as np
from scipy.special import gamma, binom
from multiFoxH import parseArgsToMultiH, firstUserMultiH
import math

pi = math.pi

def PEP_analit(user: int, L: int, params: list, theta: float, gammaBar: np.array) -> list[float]:
    points = len(gammaBar)
    # print(points)

    alpha = params[0]
    mu = params[1]
    ms = params[2]
    z = params[3]

    Psi = (mu/(ms-1)) ** (1/alpha)
    Xi = (Psi*( (z*(2**(1/2)))**alpha )) / (((theta**2) * gammaBar * (z**2 + 2))**(alpha/2))

    if user == 0:
        # first user
        sum = 0
        preH = (alpha*L) / (4*(pi**(1/2)))
        for k in range(L-1):
            binomial = binom(L-1, k)
            signal = (-1)**k
            frac = ((z**2) / (alpha*gamma(mu)*gamma(ms)))**(k+1)

            H = firstUserMultiH(params, k, Xi, points)

            sum = sum + (binomial*signal*frac*H)
        
        PEP = preH * sum
    elif user != L:
        # l_th user
        sum = 0
        gamma_ratio = (gamma(L+1)) / (gamma(user)*gamma(L-user+1))
        preH = (alpha*gamma_ratio) / (4*(pi**(1/2)))
        for k in range(L-1):
            binomial = binom(L-user, k)
            signal = (-1)**k
            frac = ((z**2) / (alpha*gamma(mu)*gamma(ms)))**(user+k)

            H = parseArgsToMultiH(params, user+k, Xi, points)

            sum = sum + (binomial*signal*frac*H)
        
        PEP = preH * sum
    else:
        # last user
        sum = 0
        frac = ((z**2) / (gamma(mu)*gamma(ms)))**L
        preH = ((alpha**(1-L))*L*frac) / (4*(pi**(1/2)))
        
        H = parseArgsToMultiH(params, L, Xi, points)
        
        PEP = preH * H

    return PEP