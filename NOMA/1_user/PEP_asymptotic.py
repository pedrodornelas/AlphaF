import numpy as np
from scipy.special import gamma, beta
import math

pi = math.pi

def PEP_asymptotic(L: int, params: list[float], theta: float, gammaBar: np.array) -> list[float]:    
    alpha = params[0]
    mu = params[1]
    ms = params[2]
    z = params[3]

    if mu < ((z**2) / alpha):
        frac1 = gamma(1/2 + (alpha*mu*L)/2) / (2 * pi**(1/2))
        frac2 = ((gamma(mu+ms) * z**2) / (alpha*mu*((z**2 / alpha)-mu)*gamma(mu)*gamma(ms)))**L
        Psi = (mu/(ms-1)) ** (1/alpha)
        frac3 = ((Psi * (2**(1/2) * z)**alpha) / (((theta**2)*gammaBar*(z**2 + 2))**(alpha/2)) )**(L*mu)
        PEP = frac1*frac2*frac3
    else:
        frac1 = gamma(1/2 + (z**2 * L)/2) / (2 * pi**(1/2))
        frac2 = ((gamma(mu - (z**2 / alpha)) * gamma(ms + (z**2 / alpha))) / (gamma(mu)*gamma(ms)))**L
        Psi = (mu/(ms-1)) ** (1/alpha)
        frac3 = ((Psi * (2**(1/2) * z)**alpha) / (((theta**2) * gammaBar * (z**2 + 2))**(alpha/2)) )**((L * z**2) / alpha)
        PEP = frac1*frac2*frac3

    return PEP