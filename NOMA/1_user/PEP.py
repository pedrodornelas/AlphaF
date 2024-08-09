import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat, loadmat
import cmath

from PEP_analit import PEP_analit
from PEP_asymptotic import PEP_asymptotic

fig, ax = plt.subplots()

# L = [1,2,3,4] # number of users
L = [3] # number of users
# alpha = [1.5, 2.3]
alpha = 2
# ms = [3, 4]
ms = 50
# mu = [1.5, 1.7]
mu = 1
# z = [0.7, 7]
z = [8]

# constelation for PEP
s1 = -0.7071 + 0.7071j
s2 = -0.7071 - 0.7071j
s3 = +0.7071 + 0.7071j
s4 = +0.7071 - 0.7071j

symbols = [[s1, s2],
           [s2, s3],
           [s3, s1]]

# power alocation to users
Beta = [0.7, 0.2, 0.1]
# power transmission
P = 1

# symbols distance
delta = symbols[0][0] - symbols[0][1]
# print(delta)

sigma = 1

# for idx, a in enumerate(alpha):
#     if ms[idx] <= (4/a):
#         print("ms > 4/alpha not met. Exiting...")
#         exit()

if ms <= (4/alpha):
    print("ms > 4/alpha not met. Exiting...")
    exit()

points = 100
l_bound_dB = 0
u_bound_dB = 30
gamma_bar_dB = np.linspace(l_bound_dB, u_bound_dB, points)
gamma_bar = 10 ** (gamma_bar_dB / 10)

PEP = np.zeros((points, len(L), len(z)))
PEP_asy = np.zeros((points, len(L), len(z)))

color = ['b','g','r','m']

for i in range(len(L)):
    for k in range(len(z)):
        analit_params = [alpha, mu, ms, z[k]]
        # print(analit_params)

        prod1 = ((Beta[0]*P)**(1/2)) * ((abs(delta))**2)
        print(prod1)
        prod2 = delta*((Beta[1]*P)**(1/2) * np.conjugate(symbols[1][0]) + (Beta[2]*P)**(1/2) * np.conjugate(symbols[2][0]))
        print(prod2)

        theta = (prod1 + 2*(prod2.real)) / (2**(1/2) * sigma * abs(delta))
        # print(theta)

        PEP[:, i, k] = PEP_analit(L[i], analit_params, theta, gamma_bar)
        # print(PEP[:, i, k])
        PEP_asy[:, i, k] = PEP_asymptotic(L[i], analit_params, theta, gamma_bar)

        if k == 0:
            ax.semilogy(gamma_bar_dB, PEP[:, i, k], color=color[i], linestyle='-', label=(f'L = {L[i]}'))
        else:
            ax.semilogy(gamma_bar_dB, PEP[:, i, k], color=color[i], linestyle='-')
        
        if i == len(L)-1 and k == len(z)-1:
            ax.semilogy(gamma_bar_dB, PEP_asy[:, i, k], color='k', linestyle='--', label='Asymptotic')
        else:
            ax.semilogy(gamma_bar_dB, PEP_asy[:, i, k], color='k', linestyle='--')

    print('L = ' + str(L[i]))

ax.legend(loc='lower left')
# ax.set_ylabel("PEP 1^st User")
# ax.set_xlabel("SNR (dB)")
ax.grid(linestyle='--', axis='both', linewidth=0.5)
plt.ylim([10**(-8), 10**0])
plt.xlim([min(gamma_bar_dB), max(gamma_bar_dB)])
plt.show()

# save = str(input('Do you want save points? (Y/N): '))
# if (save == "Y" | save == "y"):
#     filename = "PEP.mat"
#     savemat(filename, dict(L=L,
#                         gamma_bar_dB=gamma_bar_dB,
#                         PEP=PEP,
#                         PEP_asy=PEP_asy,
#                         alpha=alpha,
#                         mu=mu,
#                         ms=ms,
#                         z=z))
#     print("points saved...")