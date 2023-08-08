import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.special import comb

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams["figure.figsize"] = [5,3]


# Import Files
init_state = [ "Neel" , "HalfNeel" ]
Larray = np.arange(4, 13, 2)
T = 1
J = [0.0001, 0.001, 0.01, 0.1, 1.0]
V = 3
hz = 16
alpha = [0.50, 3.00]
kick = 0.01 #0.01
num_J = len(J)
num_alpha = len(alpha)
n_points = num_J * num_alpha
numL = len(Larray)
iter2 = 20480
steps = 10000

m = 0

filename = np.empty((2, numL, num_J, num_alpha), dtype=object)
files = np.empty((2, numL, num_J, num_alpha), dtype=object)
data = np.empty((2, numL, num_J, num_alpha), dtype=object)

sigmaz = np.empty((2, numL, num_J, num_alpha), dtype=object)
Z = np.empty((2, numL, num_J, num_alpha), dtype=object)


for s in range(2):
    for i in range(num_J):
        for j in range(num_alpha):
            for q in range(numL):
                L = Larray[q]
                print([J[i], alpha[j], L])
                n_iter = iter2 * 2 ** (1 - L / 2)
                #dim_Sz0 = comb(L, L / 2)
    
    #            filename[q, i, j] = 'data/eigen/quasienergies_Swap_LR_Sz0_nspin%d_period%.2f_n_disorder%d_Jxy%.5f_Vzz%.2f_hz%.2f_alpha%.3f_kick%.2f.txt' % (
    #                L, T, n_iter, J[i], V, hz, alpha[j], kick)
                filename[s, q, i, j] = "data/dynamics/sigmaz_Swap_LR_%s_nspin%d_period%.2f_n_disorder%d_Jxy%.5f_Vzz%.2f_hz%.2f_kick%.3f_alpha%.2f.txt" % (
                    init_state[s], L, T, n_iter, J[i], V, hz, kick, alpha[j])
    
                print(filename[s, q, i, j])
                files[s, q, i, j] = np.genfromtxt(filename[s, q, i, j], skip_header=8)
                print(files[s, q, i, j])
                data[s, q, i, j] = files[s, q, i, j]
                print(data[s, q, i, j][0, :])

for s in range(2):
    for i in range(num_J):
        for j in range(num_alpha):
            for q in range(numL):
                L = Larray[q]
                print(init_state[s] ,[J[i], alpha[j], L])
                n_iter = iter2 * 2 ** (1 - L / 2)
                #dim_Sz0 = comb(L, L / 2)
                sigmaz[s, q, i, j] = data[s, q, i, j][:, 1:L+1]
                Z[s, q, i, j] = (np.sign(sigmaz[s, q, i,j][0,0::2] - sigmaz[s,q,i,j][0,1::2]) * (sigmaz[s,q,i,j][:,0::2] - sigmaz[s,q,i,j][:,1::2])).sum( axis=1 ) / L
                print(Z[s, q,i,j][0:4],"\n")

X = np.empty((2, numL, num_J, num_alpha), dtype=object)
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
        #fig, ax = plt.subplots()
        for q in range(numL):
            for s in range(2):
                L = Larray[q]
                print(init_state[s], [J[i], alpha[j], L])
    
                X[s, q, i, j] = signs * Z[s, q,i,j] / Z[s, q,i,j][0]
                string = '$L = %d$' % L
                if (s == 1) :
                    plt.plot(range(1,steps+1), X[s, q,i,j], label=string, color=to_hex(cmL(q)))
                if (s == 0) :
                    plt.plot(range(1,steps+1), X[s, q,i,j], label=string, linestyle='--', color=to_hex(cmL(q)))
    
    
        plt.xlabel('$t$', fontsize=12)
        plt.ylabel('$(-1)^t \overline{Z}(t)/Z(0)$', fontsize=12)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.xscale('log')
        plt.ylim([-0.1,1.1])
        plt.legend(loc='lower center', bbox_to_anchor=(0.5, -0.45), fancybox=True, shadow=True, ncol=numL)
        #h, l = ax.get_legend_handles_labels()
        #ph = [plt.plot([],marker="", ls="")[0]]*2
        #handles = ph + h
        #labels = ["HalfNeel:", "Neel:"] + l
        #plt.legend(handles, labels, ncol=numL+1)
        #handles = ph[:1] + h[::2] + ph[1:] + h[1::2]
        #labels = ["Title 1:"] + l[::2] + ["Title 2:"] + l[1::2]
        #leg = plt.legend(handles, labels, ncol=numL+1)
        
        #for hpack in leg._legend_handle_box.get_children():
        #    for vpack in hpack.get_children()[:1]:
        #        hpack.get_children()[0].set_width(0)


        txt = 'figures/Z_avg_Dynamics_Comparison_wrt_L_J%.5f_alpha%.2f_kick%.3f.pdf' % (J[i], alpha[j], kick)
        print(txt)
        plt.savefig(txt, dpi=600, bbox_inches='tight')
        plt.close()
        print("\n")

