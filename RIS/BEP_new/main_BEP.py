import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat, loadmat

from BEP_analit import BEP_analit
from BEP_asymptotic import BEP_asymptotic

fig, ax = plt.subplots()

L = [1,2,3,4]
N = 2
alpha = [1.5, 2.3]
ms = [4, 5]
mu = [2.5, 3.5]
# z = [[0.7, 0.72],
#      [  7,   8]]

z = [[0.85, 0.87],
     [  7,   8]]

p = 1    # BPSK modulation

for idx, a in enumerate(alpha):
    if ms[idx] <= (4/a):
        print("ms > 4/alpha not met. Exiting...")
        exit()

points = 100
l_bound_dB = 0
u_bound_dB = 30
gamma_bar_dB = np.linspace(l_bound_dB, u_bound_dB, points)
gamma_bar = 10 ** (gamma_bar_dB / 10)

analit_gamma_bar_c = (10 ** (1 / 10)) * np.ones((points, N))
analit_gamma_bar_c[:,0] = gamma_bar
# print(analit_gamma_bar_c)

BEP = np.zeros((points, len(L), len(z)))
BEP_asy = np.zeros((points, len(L), len(z)))

color = ['b','g','r','m']

for i in range(len(L)):
    for k in range(len(z)):
        analit_params = [[alpha[0], mu[0], ms[0], z[k][0]],
                         [alpha[1], mu[1], ms[1], z[k][1]]]
        # print(analit_params)

        BEP[:, i, k] = BEP_analit(L[i], N, analit_params, p, analit_gamma_bar_c)
        # print(BEP[:, i, k])
        BEP_asy[:, i, k] = BEP_asymptotic(L[i], N, analit_params, p, analit_gamma_bar_c)

        if k == 0:
            ax.semilogy(gamma_bar_dB, BEP[:, i, k], color=color[i], linestyle='-', label=(f'L = {L[i]}'))
        else:
            ax.semilogy(gamma_bar_dB, BEP[:, i, k], color=color[i], linestyle='-')
        
        if i == len(L)-1 and k == len(z)-1:
            ax.semilogy(gamma_bar_dB, BEP_asy[:, i, k], color='k', linestyle='--', label='Asymptotic')
        else:
            ax.semilogy(gamma_bar_dB, BEP_asy[:, i, k], color='k', linestyle='--')

    print(L[i])

filename = "BEP.mat"
savemat(filename, dict(L=L,
                       N=N,
                       gamma_bar_dB=gamma_bar_dB,
                       BEP=BEP,
                       BEP_asy=BEP_asy,
                       alpha=alpha,
                       mu=mu,
                       ms=ms,
                       z=z,
                       p=p))
print("points saved...")


ax.legend(loc='lower left')
ax.set_ylabel("BEP")
ax.set_xlabel("SNR (dB)")
ax.grid(linestyle='--', axis='both', linewidth=0.5)
plt.ylim([10**(-5), 10**0])
plt.xlim([min(gamma_bar_dB), max(gamma_bar_dB)])
plt.show()