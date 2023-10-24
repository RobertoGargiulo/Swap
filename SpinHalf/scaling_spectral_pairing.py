import numpy as np
import matplotlib.pyplot as plt

from numpy import log, log10, e

import csv

# Set formatting options
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['figure.figsize'] = [5,3]

import numpy as np

import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.special import comb

###### Spectral pairing for generic parameters, for different alpha ######



# Import Files
Larray = np.arange(4, 13, 2)
alpha = [0.5, 3.0]
num_alpha = len(alpha)
nparams = 11
numL = len(Larray)
m = 0



filename = 'Swap_LR_spectrum_non_zero_J_kick.txt'

data = np.genfromtxt(filename)
print(data, "\n")

#L, T, J, V, hz, kick (, alpha x 2)
params = np.empty( (num_alpha, nparams), dtype=object)
data_list = np.empty( (num_alpha,nparams),  dtype=object)
delta_pi = np.empty( (num_alpha, nparams), dtype=object)
delta_pi_err = np.empty( (num_alpha, nparams), dtype=object)

#data_list

for j in range(num_alpha):
    for i in range(nparams):
        m1 = (i+j*nparams)*numL
        m2 = m1 + numL
        print(m1, m2, j, i)
        data_list[j][i] = data[m1:m2,:]
        params[j][i] = data_list[j][i][0,1:7]
        print(j, i, "\n", data_list[j][i][:,0],"\n", params[j][i], "\n\n")
        delta_pi[j][i] = data_list[j][i][:,11] / log(10.)
        delta_pi_err[j][i] = data_list[j][i][:,12] / log(10.)

print(np.c_[delta_pi[0][0], delta_pi_err[0][0]])
print(params[1][0][[0, 3, 4]])

X = np.empty((num_alpha, nparams), dtype=object)

for j in range(num_alpha):
    m += 1
    plt.figure(m)
    for i in [0, 1, 3, 4, 9, 10]:

        J = params[j][i][0]
        kick = params[j][i][3]
        alpha = params[j][i][4]
        print([J, kick, alpha])

        X[j][i] = delta_pi[j][i]
        string = '$J = %.g$, $\epsilon = %.g$' % ( J, kick )
        plt.errorbar(Larray, delta_pi[j][i], delta_pi_err[j][i], label=string, linewidth=1, marker='o')

    plt.xticks(Larray)
    plt.xlabel('$L$', fontsize=13)
    plt.ylabel('$\ell$', fontsize=13)
    plt.legend(fontsize=9, loc='lower left', bbox_to_anchor=(0,0))
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    txt = 'figures/Log_Delta_pi_Plot_non_zero_J_kick_alpha%.2f_L%d.pdf' % (alpha, Larray[numL-1])
    print(txt)
    plt.savefig(txt, dpi=600, bbox_inches='tight')
    print("\n")

