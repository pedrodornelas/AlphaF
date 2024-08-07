import numpy as np
from scipy.special import gamma, binom
import math

pi = math.pi

def PEP_asymptotic(L: int, l: int, params: list[float], gammaBar: np.array) -> list[float]:    
    alpha = params[0]
    mu = params[1]
    ms = params[2]
    z = params[3]

    theta = 1

    if mu < ((z**2) / alpha):
        preSum = (alpha*mu*gamma(L+1)) / (4*(pi**(1/2))*gamma(l)*gamma(L-l+1))
        sum = 0

        for k in range(L-l):
            binomial = binom((L-l), k)
            signal = (-1)**k
            frac1 = gamma(1/2 + (alpha*mu/2)*(l+k)) / (l+k)
            Psi = (mu/(ms-1)) ** (1/alpha)
            frac2 = (gamma(mu+ms)*(z**2)*(Psi**mu)) / (alpha*mu*gamma(mu)*gamma(ms)*((z**2 / alpha)-mu))
            frac3 = ((z * 2**(1/2)) / (((theta**2) * gammaBar * (z**2 + 2))**(1/2)))**(alpha*mu)
            sum = sum + (binomial*signal*frac1*((frac2*frac3)**(l+k)))

        PEP = preSum * sum

    else:
        preSum = (alpha*mu*gamma(L+1)) / (4*(pi**(1/2))*gamma(l)*gamma(L-l+1))
        sum = 0

        for k in range(L-l):
            binomial = binom((L-l), k)
            signal = (-1)**k
            frac1 = gamma(1/2 + (z**2 / 2)*(l+k)) / (l+k)
            Psi = (mu/(ms-1)) ** (1/alpha)
            frac2 = (gamma(mu-(z**2 / alpha))*gamma(ms+(z**2 / alpha))*(Psi**(z**2 / alpha))) / (gamma(mu)*gamma(ms))
            frac3 = ((z * 2**(1/2)) / (((theta**2) * gammaBar * (z**2 + 2))**(1/2)))**(z**2)
            sum = sum + (binomial*signal*frac1*((frac2*frac3)**(l+k)))

        PEP = preSum * sum

    return PEP