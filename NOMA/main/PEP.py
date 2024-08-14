import matplotlib.pyplot as plt
import numpy as np
from scipy.io import savemat, loadmat
import cmath

from PEP_analit import PEP_analit
from PEP_asymptotic import PEP_asymptotic

def createTheta(user: int, L: int, Beta: list[float], P: float, sigma: float, symbols) -> float:
    sym_tx = symbols[user][0]
    sym_rx = symbols[user][1]
    # symbols distance
    delta_user = sym_tx - sym_rx

    denom = (2**(1/2)) * sigma * abs(delta_user)

    prod1 = ((Beta[user]*P)**(1/2)) * ((abs(delta_user))**2)
    # print(prod1)

    sum1 = 0
    if user != 0:
        for i in range(0, user):
            print('i = ' + str(i))
            delta_l = symbols[i][0] - symbols[i][1]
            sum1 = sum1 + ((Beta[i]*P)**(1/2))*(np.conjugate(delta_l))

    sum2 = 0
    if user != (L-1):
        for j in range(user+1, L):
            print('j = ' + str(j))
            sum2 = sum2 + (Beta[j]*P)**(1/2)*(np.conjugate(symbols[j][0]))

    sum = sum1 + sum2
    prod2 = delta_user*sum
    # print(prod2)

    final_prod = (prod1 + 2*(prod2.real)) / denom
    return final_prod

fig, ax = plt.subplots()

# L = [1,2,3,4] # number of users
L = 3 # number of users
# alpha = [1.5, 2.3]
alpha = 2
# ms = [3, 4]
ms = 50
# mu = [1.5, 1.7]
mu = 1
# z = [0.7, 7]
z = [8]

# constelation for PEP
# s1 = (-0.7071 + 0.7071j)
# s2 = (-0.7071 - 0.7071j)
# s3 = (+0.7071 + 0.7071j)
# s4 = (+0.7071 - 0.7071j)

s1 = (-1 - 1j)
s2 = (-1 + 1j)
s3 = (+1 - 1j)
s4 = (+1 + 1j)

symbols = [[s2, s1],
           [s2, s1],
           [s2, s1]]

print(symbols)
# symbols = np.multiply(L**(1/2), symbols)
# print(symbols)

# print(np.multiply((L**(1/2)), symbols))

# power alocation to users
Beta = [0.7, 0.2, 0.1]
# power transmission
P = 1

sigma = 1

# for idx, a in enumerate(alpha):
#     if ms[idx] <= (4/a):
#         print("ms > 4/alpha not met. Exiting...")
#         exit()

if ms <= (4/alpha):
    print("ms > 4/alpha not met. Exiting...")
    exit()

points = 100
l_bound_dB = -30
u_bound_dB = 20
gamma_bar_dB = np.linspace(l_bound_dB, u_bound_dB, points)
gamma_bar = 10 ** (gamma_bar_dB / 10)
# gamma_bar = (L**(1/2)) * (10 ** (gamma_bar_dB / 10))


PEP = np.zeros((points, L, len(z)))
PEP_asy = np.zeros((points, L, len(z)))

color = ['b','g','r','m']

for l in range(L):
    print('user = ' + str(l))
    for k in range(len(z)):
        analit_params = [alpha, mu, ms, z[k]]
        # print(analit_params)

        theta = createTheta(l, L, Beta, P, sigma, symbols)
        # print(theta)

        PEP[:, l, k] = PEP_analit(l, L, analit_params, theta, gamma_bar)
        # print(PEP[:, l, k])
        # PEP_asy[:, l, k] = PEP_asymptotic(l, L, analit_params, theta, gamma_bar)

        if k == 0:
            ax.semilogy(gamma_bar_dB, PEP[:, l, k], color=color[l], linestyle='-', label=(f'U = {l}'))
        else:
            ax.semilogy(gamma_bar_dB, PEP[:, l, k], color=color[l], linestyle='-')
        
        # if l == L-1 and k == len(z)-1:
        #     ax.semilogy(gamma_bar_dB, PEP_asy[:, l, k], color='k', linestyle='--', label='Asymptotic')
        # else:
        #     ax.semilogy(gamma_bar_dB, PEP_asy[:, l, k], color='k', linestyle='--')

ax.legend(loc='lower left')
# ax.set_ylabel("PEP 1^st User")
# ax.set_xlabel("SNR (dB)")
ax.grid(linestyle='--', axis='both', linewidth=0.5)
plt.ylim([10**(-8), 10**0])
plt.xlim([min(gamma_bar_dB), max(gamma_bar_dB)])
plt.show()

# save = str(input('Do you want save points? (Y/N): '))
# if (save == "Y" or save == "y"):
#     filename = "PEP.mat"
#     savemat(filename, dict(L=L,
#                         gamma_bar_dB=gamma_bar_dB,
#                         PEP=PEP,
#                         PEP_asy=PEP_asy,
#                         alpha=alpha,
#                         mu=mu,
#                         ms=ms,
#                         z=z))
#     print("saved points...")
# else:
#     print("discarded points")