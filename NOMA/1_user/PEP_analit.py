import numpy as np
from scipy.special import gamma
from multiFoxH import parseArgsToMultiH
import math

pi = math.pi

def PEP_analit(L: int, params: list, theta: float, gammaBar: np.array) -> list[float]:
    points = len(gammaBar)
    # print(points)

    alpha = params[0]
    mu = params[1]
    ms = params[2]
    z = params[3]

    preH = ( (L * ((alpha)**(1-L))) / (4*( (pi)**(1/2) ))) * (( (z**2) / (gamma(mu)*gamma(ms)))**L)

    Psi = (mu/(ms-1)) ** (1/alpha)
    Xi = (Psi*( (z*(2**(1/2)))**alpha )) / (((theta**2) * gammaBar * (z**2 + 2))**(alpha/2))
    # print(Xi)

    H = parseArgsToMultiH(params, L, Xi, points)
    # print(len(H))
    # print(H)

    PEP = preH * H
    # print(len(PEP))
    # print(PEP)

    return PEP