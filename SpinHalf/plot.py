import numpy as np

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.special import comb

# Import Files
Larray = np.arange(4, 15, 2)
T = 1
J = [0.00001, 0.0001, 0.001, 0.01, 0.1]
V = 3
hz = 16
kick = [0.00]
alpha = 0.50
num_J = len(J)
num_kick = len(kick)
n_points = num_J * num_kick
numL = len(Larray)
iter2 = 20480

m = 0

filename = np.empty((numL, num_J, num_kick), dtype=object)
files = np.empty((numL, num_J, num_kick), dtype=object)
data = np.empty((numL, num_J, num_kick), dtype=object)
arr = np.empty((numL, 1), dtype=object)

QE = np.empty((numL, num_J, num_kick), dtype=object)
E_MBL = np.empty((numL, num_J, num_kick), dtype=object)

edges = np.empty((numL, num_J, num_kick), dtype=object)
counts = np.empty((numL, num_J, num_kick), dtype=object)

for i in range(num_J):
    for j in range(num_kick):
        for q in range(numL):
            L = Larray[q]
            print([J[i], kick[j], L])
            n_iter = iter2 * 2 ** (1 - L / 2)
            dim_Sz0 = comb(L, L / 2)

            filename[q, i, j] = 'data/eigen/quasienergies_Swap_LR_Sz0_nspin%d_period%.2f_n_disorder%d_Jxy%.5f_Vzz%.2f_hz%.2f_kick%.3f_alpha%.2f.txt' % (
                L, T, n_iter, J[i], V, hz, kick[j], alpha)

            print(filename[q, i, j])
            files[q, i, j] = np.genfromtxt(filename[q, i, j], skip_header=7)
            print(files[q, i, j])
            data[q, i, j] = files[q, i, j]
            print(data[q, i, j][0, :])

for i in range(num_J):
    for j in range(num_kick):
        for q in range(numL):
            L = Larray[q]
            print([J[i], kick[j], L])
            n_iter = iter2 * 2 ** (1 - L / 2)
            dim_Sz0 = comb(L, L / 2)
            QE[q, i, j] = data[q, i, j][:, 2]
            E_MBL[q, i, j] = data[q, i, j][:, 3]

from numpy import pi

def pi_pair(qe):
    pi_pair = np.zeros(np.shape(qe), dtype=int)
    for i in range(len(qe)):
        val = (( qe[i] + 2*pi) % (2*pi)) - pi
        #print(qe[i], (qe[i]) % (2*pi) - pi, val)
        pi_pair[i] = np.argmin( np.abs( np.exp( 1j*(qe - val) ) - 1 ) )

    return pi_pair

def pi_pairs(qe):
    pi_pair = np.zeros(np.shape(qe), dtype=int)
    for i in range(len(qe)):
        val = (( qe[i] + 2*pi) % (2*pi)) - pi
        #print(qe[i], (qe[i]) % (2*pi) - pi, val)
        pi_pair[i] = np.argmin( np.abs( np.exp( 1j*(qe - val) ) - 1 ) )
        #print(np.abs( np.exp( 1j*(qe - val) ) - 1 ))

        #print(qe[pi_pair])
        #print(qe)
        pair = np.abs(np.abs(qe[pi_pair] - qe) - pi)
    return pair
    

#print(np.c_[qe[:,None], qe[pi_pair,None], np.abs(np.abs(qe[pi_pair,None] - qe[:,None]) - pi), \
#           np.abs(np.abs(qe[(pi_pair + 1)%len(qe),None] - qe[:,None]) - pi) ] )

qe = 2*pi * np.random.random( (1000,) ) - pi
qe = np.sort(qe)
plt.plot(qe)
plt.plot(qe[pi_pair(qe)])

txt = 'figures/figure2.png' #% (J[i], kick[j])
print(txt)
plt.savefig(txt, dpi=600)
#plt.close()
print("\n")
