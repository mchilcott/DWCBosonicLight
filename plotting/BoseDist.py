import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve 


def N_BE(mu, omega):
    # mu in 1/kT
    # omega in hbar/kT

    n = np.arange(0, 1000)
    return (n+1) / (np.exp(n*omega - mu) - 1)

def find_mu (mu, omega, N):
    return np.sum(N_BE(mu, omega)) - N


N = 2000
Tem = np.linspace(1, 70, 1000) # in units k/(hbar * omega)
n_0 = np.zeros(Tem.size)

for i in range(0, Tem.size):

    T = Tem[i];
    mu = fsolve(func=find_mu, x0=-0.001, args=(1/T, N))
    n_0[i] = N_BE(mu, 1/T)[0]




T_C = np.sqrt(6*N)/np.pi


fig, ax1 = plt.subplots()
ax1.plot(Tem/T_C, n_0/N, label="Numeric")
ax1.set_xlabel("Temperature [$T/T_C$]")
ax1.set_ylabel("Ground state (BEC) fraction")
ax1.axvline(1, c="tab:orange")

ax2 = ax1.twiny()
ax2.set_xlabel(r"Temperature $\left[\frac{k_B}{\hbar\omega}\right]$,  N=2000")
ax2.set_xlim(np.stack(ax1.get_xlim()) * T_C)

plt.tight_layout()
plt.savefig("../presentation/gfx/BECFrac.svg")
plt.show()
