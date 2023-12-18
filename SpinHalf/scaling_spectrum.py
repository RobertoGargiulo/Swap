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
Larray = np.arange(4, 15, 2)
alpha = [0.5, 3.0]
num_alpha = len(alpha)
numL = len(Larray)
m = 0
J = [0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0]
num_J = len(J)

filename = "sort_col2_Swap_LR_spectrum_kick0.txt"

data = np.genfromtxt(filename)
print(data, "\n")

#L, T, J, V, hz, kick (, alpha x 2)
params = np.empty( (num_alpha, numL), dtype=object)
data_list = np.empty( (num_alpha,numL),  dtype=object)
delta_pi = np.empty( (num_alpha, numL), dtype=object)
delta_pi_err = np.empty( (num_alpha, numL), dtype=object)
gap = np.empty( (num_alpha, numL), dtype=object)
gap_err = np.empty( (num_alpha, numL), dtype=object)

#data_list

for j in range(num_alpha):
    for i in range(numL):
        m1 = (i+j*numL)*num_J
        m2 = m1 + num_J
        print(m1, m2, j, i)
        data_list[j][i] = data[m1:m2,:]
        params[j][i] = data_list[j][i][0,0:7]
        print(j, i, "\n", data_list[j][i][:,1],"\n", params[j][i], "\n\n")
        delta_pi[j][i] = data_list[j][i][:,11] / log(10.)
        delta_pi_err[j][i] = data_list[j][i][:,12] / log(10.)
        gap[j][i] = data_list[j][i][:,7]
        gap_err[j][i] = data_list[j][i][:,8]

print(np.c_[delta_pi[0][0], delta_pi_err[0][0]])
print(np.c_[gap[0][0], gap_err[0][0]])
print(params[0][0][[0,1,2,5]])


for j in range(num_alpha):
    m += 1
    plt.figure(m)
    for i in range(numL):

        L = Larray[i]
        kick = params[j][i][4]
        alpha = params[j][i][5]
        print([L, kick, alpha])

        string = '$L = %d$' % ( L )
        plt.errorbar(J, delta_pi[j][i], delta_pi_err[j][i], label=string, linewidth=1, marker='o')

    #plt.xticks(J)
    plt.tick_params(labelbottom=False)
    plt.xscale('log')
    #plt.xlabel('$J$', fontsize=12)
    plt.ylabel('$\ell_\Delta$', fontsize=12)
    plt.legend(loc='lower right', fancybox=True, shadow=True, fontsize = 10)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    txt = 'figures/Log_Delta_pi_Plot_kick%.2f_alpha%.2f_L%d.pdf' % (kick, alpha, Larray[numL-1])
    print(txt)
    plt.savefig(txt, dpi=600, bbox_inches='tight')
    print("\n")

for j in range(num_alpha):
    m += 1
    plt.figure(m)
    for i in range(numL):

        L = Larray[i]
        kick = params[j][i][4]
        alpha = params[j][i][5]
        print([L, kick, alpha])

        string = '$L = %d$' % ( L )
        plt.errorbar(J, gap[j][i], gap_err[j][i], label=string, linewidth=1, marker='o')

    plt.axhline(y=0.386, color='black', linestyle='--', label='Poisson')
    plt.axhline(y=0.5269, color='black', linestyle=':', label='COE')
    plt.tick_params(labelbottom=False)
    #plt.xticks(J)
    plt.xticks([])
    plt.xscale('log')
   # plt.xlabel('$J$', fontsize=12)
    plt.ylabel('$\langle r \\rangle$', fontsize=12)
    plt.ylim([0.33,0.6])
    plt.legend(loc='upper center', fancybox=True, shadow=True, fontsize = 10, ncol=numL-2)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    txt = 'figures/Gap_Ratio_Plot_kick%.2f_alpha%.2f_L%d.pdf' % (kick, alpha, Larray[numL-1])
    print(txt)
    plt.savefig(txt, dpi=600, bbox_inches='tight')
    print("\n")
