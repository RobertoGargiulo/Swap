import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.special import comb
from matplotlib import gridspec

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams["figure.figsize"] = [5,6]


# Import Files
init_state = "UpZero"
Larray = np.arange(4, 11, 2)
T = 1
J = [0.0001, 0.001, 0.01, 0.1, 1.0]
V = 3
hz = 16
alpha = [0.50, 3.00]
kick = 0.00 #0.01
num_J = len(J)
num_alpha = len(alpha)
n_points = num_J * num_alpha
numL = len(Larray)
iter2 = 20480
steps = 10000

m = 0

filename = np.empty(( numL, num_J, num_alpha), dtype=object)
files = np.empty(( numL, num_J, num_alpha), dtype=object)
data = np.empty(( numL, num_J, num_alpha), dtype=object)

sigmaz = np.empty(( numL, num_J, num_alpha), dtype=object)
Z = np.empty(( numL, num_J, num_alpha), dtype=object)


for i in range(num_J):
    for j in range(num_alpha):
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])
            n_iter = iter2 * 2 ** (1 - L / 2)
            #dim_Sz0 = comb(L, L / 2)

#            filename[q, i, j] = 'data/eigen/quasienergies_Swap_LR_Sz0_nspin%d_period%.2f_n_disorder%d_Jxy%.5f_Vzz%.2f_hz%.2f_alpha%.3f_kick%.2f.txt' % (
#                L, T, n_iter, J[i], V, hz, alpha[j], kick)
            filename[ q, i, j] = "data/dynamics/sigmaz_Swap_LR_%s_nspin%d_steps%d_period%.2f_n_disorder%d_J%.5f_V%.2f_hz%.2f_kick%.3f_alpha%.2f.txt" % (
                init_state, L, steps, T, n_iter, J[i], V, hz, kick, alpha[j])

            print(filename[ q, i, j])
            files[ q, i, j] = np.genfromtxt(filename[ q, i, j], skip_header=8)
            print(files[ q, i, j])
            data[ q, i, j] = files[ q, i, j]
            print(data[ q, i, j][0, :])

for i in range(num_J):
    for j in range(num_alpha):
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])
            n_iter = iter2 * 2 ** (1 - L / 2)
            #dim_Sz0 = comb(L, L / 2)
            sigmaz[ q, i, j] = data[ q, i, j][:, 1:L+1]
            Z[ q, i, j] = (np.sign(sigmaz[ q, i,j][0,0::2] - sigmaz[q,i,j][0,1::2]) * (sigmaz[q,i,j][:,0::2] - sigmaz[q,i,j][:,1::2])).sum( axis=1 ) / L
            print(Z[ q,i,j][0:4],"\n")

X = np.empty((numL, num_J, num_alpha), dtype=object)
Y = np.empty((numL, num_J, num_alpha), dtype=object)
signs = np.empty((steps,))
signs[::2] = +1
signs[1::2] = -1

plt.rcParams["axes.prop_cycle"] = plt.cycler("color", plt.cm.rainbow(np.linspace(0,1,numL)))

import matplotlib.cm as cm
from matplotlib.colors import to_hex

cmL = cm.get_cmap('turbo', numL)

for i in range(num_J):
    for j in range(num_alpha):
        m += 1
        plt.figure(m)
        fig, (ax1, ax2) = plt.subplots(nrows=2, sharex=True)
        plt.sca(ax2)
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])

            X[q, i, j] = signs * (sigmaz[q,i,j][:,2*(L//4)-1] - 0.5) + 0.5
            Y[q, i, j] = signs * (sigmaz[q,i,j][:,2*(L//4)] - 0.5) + 0.5
            print(sigmaz[q,i,j][0,:])
            print(sigmaz[q,i,j][0,2*(L//4)], sigmaz[q,i,j][0,L-1])
            string = '$L = %d$' % L
            plt.plot(range(1,steps+1), X[q,i,j], label=string, color=to_hex(cmL(q)))
            plt.plot(range(1,steps+1), Y[q,i,j], label=string, linewidth=0.8, linestyle='--', color=to_hex(cmL(q)))


        plt.xlabel('$t$', fontsize=12)
        plt.xticks(fontsize=10)
        plt.ylabel('$(-1)^t(\langle\overline{S_k^z}(t)\\rangle-0.5)+0.5$', fontsize=12)
        plt.yticks(fontsize=10)
        plt.xscale('log')
        plt.ylim([-0.1,1.1])
        plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.45), fancybox=True, shadow=True, ncol=numL)

        plt.sca(ax1)
        plt.subplots_adjust(hspace=.2)
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])

            X[q, i, j] = signs * Z[q,i,j] / Z[q,i,j][0]
            string = '$L = %d$' % L
            plt.plot(range(1,steps+1), X[q,i,j], label=string)


        plt.ylabel('$(-1)^t \overline{Z}(t)/Z(0)$', fontsize=12)
        plt.yticks(fontsize=10)
        plt.xscale('log')
        plt.ylim([-0.1,1.1])
        plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.2), fancybox=True, shadow=True, ncol=numL)

        txt = 'figures/sigmaz_Dynamics_two_panel_%s_wrt_L_J%.5f_alpha%.2f_kick%.3f.pdf' % (init_state, J[i], alpha[j], kick)
        print(txt)
        plt.savefig(txt, dpi=600, bbox_inches='tight')
        print("\n")

#end = 20
#for i in range(num_J):
#    for j in range(num_alpha):
#        m += 1
#        plt.figure(m)
#        for q in range(numL):
#            L = Larray[q]
#            print([J[i], alpha[j], L])
#
#            #X[q, i, j] = signs * (sigmaz[q,i,j][:,2*(L//4)] - 0.5) + 0.5
#            #Y[q, i, j] = signs * (sigmaz[q,i,j][:,L-1] - 0.5) + 0.5
#            X[q, i, j] = sigmaz[q,i,j][:,2*(L//4)]
#            Y[q, i, j] = sigmaz[q,i,j][:,L-1]
#            print(sigmaz[q,i,j][0,:])
#            print(sigmaz[q,i,j][0,2*(L//4)], sigmaz[q,i,j][0,L-1])
#            string = '$L = %d$' % L
#            plt.plot(range(1,end+1), X[q,i,j][0:end], label=string, color=to_hex(cmL(q)))
#            #plt.plot(range(1,end+1), Y[q,i,j][0:end-1], label=string, linewidth=0.8, linestyle='--', color=to_hex(cmL(q)))
#
#
#        plt.xlabel('$t$', fontsize=12)
#        plt.ylabel('$(-1)^t(\overline{\sigma_k^z}(t)-0.5)+0.5$', fontsize=12)
#        plt.xticks(fontsize=10)
#        plt.yticks(fontsize=10)
#        plt.xscale('log')
#        plt.ylim([-0.1,1.1])
#        plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.45), fancybox=True, shadow=True, ncol=numL)
#
#        txt = 'figures/sigmaz_PARTIAL_avg_Dynamics_%s_Comparison_wrt_L_J%.5f_alpha%.2f_kick%.3f.pdf' % (init_state, J[i], alpha[j], kick)
#        print(txt)
#        plt.savefig(txt, dpi=600, bbox_inches='tight')
#        print("\n")

