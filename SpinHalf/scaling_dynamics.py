import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.special import comb
from numpy import log

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams["figure.figsize"] = [5,3]


# Import Files
init_state = "Neel"
Larray = np.arange(4, 13, 2)
T = 1
J = [1.0] #[0.01, 0.1]
V = 3
hz = 16
alpha = [0.50, 3.00]
kick = 0.00 #0.01
num_J = len(J)
num_alpha = len(alpha)
n_points = num_J * num_alpha
numL = len(Larray)
iter2 = 10240
steps = 1000000
skip = 1

m = 0

filename = np.empty((numL, num_J, num_alpha), dtype=object)
files = np.empty((numL, num_J, num_alpha), dtype=object)
data = np.empty((numL, num_J, num_alpha), dtype=object)
arr = np.empty((numL, 1), dtype=object)

times = np.empty((numL, num_J, num_alpha), dtype=object)
sigmaz = np.empty((numL, num_J, num_alpha), dtype=object)
Z = np.empty((numL, num_J, num_alpha), dtype=object)
n_periods = np.empty((numL, num_J, num_alpha), dtype=int)


for i in range(num_J):
    for j in range(num_alpha):
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])
            n_iter = iter2 * 2 ** (1 - L / 2)

            if J[i] >= 0.25:
                n_periods[q,i,j] = 2**0
            elif J[i] < 0.25 and J[i] > 0.05:
                n_periods[q,i,j] = 2**( max(0, int(log(0.42 * (3.1*J[i])**(-0.38*L) ) / log(2))+1) )
            elif J[i] <= 0.05 and J[i] > 0.01:
                n_periods[q,i,j] = 2**( max(0, int(log(0.42 * (3.1*J[i])**(-0.38*L) ) / log(2))+2) )
                if L >= 8:
                    n_periods[q,i,j] = 2*n_periods[q,i,j]
            elif J[i] <= 0.01:
                n_periods[q,i,j] = 2**( max(0, int(log(0.42 * (3.1*J[i])**(-0.38*L) ) / log(2))+3) )
                if L >= 8:
                    n_periods[q,i,j] = 2*n_periods[q,i,j]

            filename[q, i, j] = "data/dynamics/sigmaz_Swap_LR_%s_nspin%d_steps%d_period%.2f_n_disorder%d_n_periods%d_Jxy%.5f_Vzz%.2f_hz%.2f_kick%.3f_alpha%.2f.txt" % (
                init_state, L, steps, T, n_iter, n_periods[q,i,j], J[i], V, hz, kick, alpha[j])

            print(filename[q, i, j])
            files[q, i, j] = np.genfromtxt(filename[q, i, j], skip_header=9)
            print(files[q, i, j])
            data[q, i, j] = files[q, i, j]
            print(data[q, i, j][0, :])

for i in range(num_J):
    for j in range(num_alpha):
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])
            n_iter = iter2 * 2 ** (1 - L / 2)
            #dim_Sz0 = comb(L, L / 2)
            times[q, i, j] = data[q, i, j][:, 0]
            sigmaz[q, i, j] = data[q, i, j][:, 1:L+1:skip]
            Z[q, i, j] = (np.sign(sigmaz[q,i,j][0,0::2] - sigmaz[q,i,j][0,1::2]) * (sigmaz[q,i,j][:,0::2] - sigmaz[q,i,j][:,1::2])).sum( axis=1 ) / L
            print(Z[q,i,j][0:4],"\n")

X = np.empty((numL, num_J, num_alpha), dtype=object)
Y = np.empty((numL, num_J, num_alpha), dtype=object)
T = np.empty((numL, num_J, num_alpha), dtype=object)


for i in range(num_J):
    for j in range(num_alpha):
        m += 1
        plt.figure(m)
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])

            T[q,i,j] = times[q,i,j]
            X[q, i, j] = (-1)**(T[q,i,j]-1) * Z[q,i,j] / Z[q,i,j][0]
            string = '$L = %d$' % L
            plt.plot(T[q,i,j], X[q,i,j], label=string)


        plt.xlabel('$t$', fontsize=12)
        plt.ylabel('$(-1)^t \overline{Z}(t)$', fontsize=12)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.xscale('log')
        plt.ylim([-0.1,1.1])
        plt.legend()

        txt = 'figures/Z_avg_Dynamics_%s_wrt_L_J%.5f_alpha%.2f_kick%.3f_long_w_skip.pdf' % (init_state, J[i], alpha[j], kick)
        print(txt)
        plt.savefig(txt, dpi=600, bbox_inches='tight')
        print("\n")

"""


import matplotlib.cm as cm
from matplotlib.colors import to_hex

cmL = cm.get_cmap('turbo', numL)

for i in range(num_J):
    for j in range(num_alpha):
        m += 1
        plt.figure(m)
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])

            X[q, i, j] =  * np.sign(sigmaz[q,i,j][0,int(L/4)]) * sigmaz[q,i,j][:,int(L/4)]
            Y[q,i,j] = sigmaz[q,i,j][:,int(L-L/4)]
            string = '$L = %d$' % L
            plt.plot(range(1,steps+1), X[q,i,j], label=string, color=to_hex(cmL(q)))
            plt.plot(range(1,steps+1), Y[q,i,j], label=string, linestyle='--', color=to_hex(cmL(q)))


        plt.xlabel('$t$', fontsize=12)
        plt.ylabel('$\langle\overline{\sigma_k^z}(t)\\rangle$', fontsize=12)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.xscale('log')
        plt.ylim([-0.1,1.1])
        plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.45), fancybox=True, shadow=True, ncol=numL)

        txt = 'figures/sigmaz_avg_Dynamics_%s_Comparison_wrt_L_J%.5f_alpha%.2f_kick%.3f_steps%d_long_w_skip.pdf' % (init_state, J[i], alpha[j], kick, steps)
        print(txt)
        plt.savefig(txt, dpi=600, bbox_inches='tight')
        print("\n")

"""
