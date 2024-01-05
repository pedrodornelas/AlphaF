import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat, loadmat

from OP_analit import OP_analit
from OP_asymptotic import OP_asymptotic

fig, ax = plt.subplots()

L = [1,2,3]
N = 2
alpha = [1.5, 2.3]
ms = [3, 4]
mu = [1.5, 1.7]
z = [[0.7, 0.8],
     [  7,   8]]

gamma_th_dB = 5
gamma_th = 10 ** (gamma_th_dB / 10)

for idx, a in enumerate(alpha):
    if ms[idx] <= (4/a):
        print("ms > 4/alpha not met. Exiting...")
        exit()

points = 100
l_bound_dB = 0
u_bound_dB = 50
gamma_bar_dB = np.linspace(l_bound_dB, u_bound_dB, points)
gamma_bar = 10 ** (gamma_bar_dB / 10)

analit_gamma_bar_c = (10 ** (1 / 10)) * np.ones((points, N))
analit_gamma_bar_c[:,0] = gamma_bar
# print(analit_gamma_bar_c)

OP = np.zeros((points, len(L), len(z)))
OP_asy = np.zeros((points, len(L), len(z)))

color = ['b','g','r','m']

for i in range(len(L)):
    for k in range(len(z)):
        analit_params = [[alpha[0], mu[0], ms[0], z[k][0]],
                         [alpha[1], mu[1], ms[1], z[k][1]]]
        # print(analit_params)

        OP[:, i, k] = OP_analit(L[i], N, analit_params, gamma_th, analit_gamma_bar_c)
        # print(OP[:, i, k])
        OP_asy[:, i, k] = OP_asymptotic(L[i], N, analit_params, gamma_th, analit_gamma_bar_c)

        if k == 0:
            ax.semilogy(gamma_bar_dB, OP[:, i, k], color=color[i], linestyle='-', label=(f'L = {L[i]}'))
        else:
            ax.semilogy(gamma_bar_dB, OP[:, i, k], color=color[i], linestyle='-')
        
        if i == len(L)-1 and k == len(z)-1:
            ax.semilogy(gamma_bar_dB, OP_asy[:, i, k], color='k', linestyle='--', label='Asymptotic')
        else:
            ax.semilogy(gamma_bar_dB, OP_asy[:, i, k], color='k', linestyle='--')

    print(L[i])

filename = "OP.mat"
savemat(filename, dict(L=L,
                       N=N,
                       gamma_bar_dB=gamma_bar_dB,
                       OP=OP,
                       OP_asy=OP_asy))
print("points saved...")


ax.legend(loc='lower left')
ax.set_ylabel("OP")
ax.set_xlabel("SNR (dB)")
ax.grid(linestyle='--', axis='both', linewidth=0.5)
plt.ylim([10**(-5), 10**0])
plt.xlim([min(gamma_bar_dB), max(gamma_bar_dB)])
plt.show()