import numpy as np
from scipy.special import gamma, beta
import math

pi = math.pi

def PEP_asymptotic(user: int, L: int, params: list[float], theta: float, gammaBar: np.array) -> list[float]:    
    alpha = params[0]
    mu = params[1]
    ms = params[2]
    z = params[3]

    real_user = user+1

    if mu < ((z**2) / alpha):
        frac1 = 1 / (2*real_user* pi**(1/2))
        frac2 = (gamma(L+1)*gamma((1 + alpha*mu*real_user)/2)) / (gamma(L)*gamma(L-real_user+1))
        Psi = (mu/(ms-1)) ** (1/alpha)
        frac3 = (Psi**mu * z**2 * gamma(mu+ms)) / (alpha*mu*gamma(mu)*gamma(ms)*(((z**2)/alpha)-mu))
        frac4 = (z * 2**(1/2)) / ((theta**2 * gammaBar * (z**2 + 2))**(1/2))
        PEP = frac1*frac2*((frac3*((frac4)**(alpha*mu)))**real_user)
    else:
        frac1 = 1 / (2*real_user* pi**(1/2))
        frac2 = (gamma(L+1)*gamma((1 + (z**2)*real_user)/2)) / (gamma(L)*gamma(L-real_user+1))
        Psi = (mu/(ms-1)) ** (1/alpha)
        frac3 = (Psi**((z**2)/alpha) * gamma(mu-((z**2)/alpha)) * gamma(ms+((z**2)/alpha))) / (gamma(mu)*gamma(ms))
        frac4 = (z * 2**(1/2)) / ((theta**2 * gammaBar * (z**2 + 2))**(1/2))
        PEP = frac1*frac2*((frac3*((frac4)**(z**2)))**real_user)

    return PEP