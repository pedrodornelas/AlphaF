import numpy as np
from scipy.special import gamma

def OP_asymptotic(L: int, N: int, params: list, gamma_th: float, gammaBar: np.array) -> list:
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
    # print(np.size(Xi))
    preH = (preH ** L) / 2
    # print(preH)
    
    alphas = np.array(params)[:, 0]
    mus = np.array(params)[:, 1]
    mss = np.array(params)[:, 2]
    zs = np.array(params)[:, 3]
    
    an = np.concatenate(([1], (1 - mss)))
    # print(an)

    An = np.concatenate(([1], (1 / alphas)))
    # print(An)

    ap = np.concatenate([((zs**2) / alphas + 1)])
    # print(ap)

    Ap = np.concatenate([(1 / alphas)])
    # print(Ap)

    bm = np.concatenate([mus, (zs**2 / alphas)])
    # print(bm)

    Bm = np.concatenate([(1 / alphas), (1 / alphas)])
    # print(Bm)

    bq = np.array([])
    Bq = np.array([])

    div = []
    for i in range(len(bm)):
        div.append(bm[i]/Bm[i])
    Ui = min(div)
    j = div.index(Ui)

    Bc = 1
    Bc = Bm[j] ** L

    gamma1 = 1
    for i in range(len(bm)):
        if i != j:
            gamma1 = gamma1*gamma(bm[i] - Ui*Bm[i])
    
    gamma2 = 1
    for i in range(len(an)):
        gamma2 = gamma2*gamma(1 - an[i] + Ui*An[i])
    
    gamma3 = 1
    for i in range(len(bq)):
        gamma3 = gamma3*gamma(1 - bq[i] + Ui*Bq[i])

    gamma4 = 1
    for i in range(len(ap)):
        gamma4 = gamma4*gamma(ap[i] - Ui*Ap[i])

    phiU = ((gamma1*gamma2)/(gamma3*gamma4)) ** L

    gamma_sum1 = gamma((Ui/2) * L)
    gamma_sum2 = gamma(Ui * L)
    gamma_sum3 = gamma(1 + ((Ui/2) * L))
    psiU = gamma_sum1 / (gamma_sum2*gamma_sum3)

    OP = (preH * psiU * phiU * (Xi**(Ui*L))) / Bc

    return OP