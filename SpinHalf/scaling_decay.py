import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from matplotlib.colors import to_hex
import matplotlib.cm as cm
from matplotlib.backends.backend_pdf import PdfPages
from numpy import sqrt, log10, log

print("-----------------------------------")

# Set formatting options
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['figure.figsize'] = [5,3]


# Import Files
Larray = np.arange(4, 13, 2)
T = 1
J = [0.05, 0.1, 0.5, 2.0] 
V = 3
hz = 16
kick = 0.00
alpha = [0.50, 3.00]
num_J = len(J)
num_alpha = len(alpha)
n_points = num_J * num_alpha
numL = len(Larray)
iter2 = 20480
steps = 1000000


m = 0

filename = np.empty((numL, num_J, num_alpha), dtype=object)
files = np.empty((numL, num_J, num_alpha), dtype=object)
data = np.empty((numL, num_J, num_alpha), dtype=object)

tau = np.empty((numL, num_J, num_alpha), dtype=object)
n_periods = np.empty((numL, num_J, num_alpha))


for i in range(num_J):
    for j in range(num_alpha):
        for q in range(numL):
            L = Larray[q]
            n_iter = iter2 * 2 ** (1 - L / 2)
            
            if J[i] >= 0.25:
                n_periods = 2**0
            elif J[i] < 0.25 and J[i] > 0.05:
                n_periods = 2**( max(0, int(log(0.42 * (3.1*J[i])**(-0.38*L) ) / log(2))+1) )
            elif J[i] <= 0.05:
                n_periods = 2**( max(0, int(log(0.42 * (3.1*J[i])**(-0.38*L) ) / log(2))+2) )
                if L >= 8:
                    n_periods = 2*n_periods
            print([J[i], alpha[j], L, n_periods])
            


            filename[q, i, j] = 'data/dynamics/decay_sigmaz_Swap_LR_Neel_nspin%d_steps%d_period%.2f_n_disorder%d_n_periods%d_Jxy%.5f_Vzz%.2f_hz%.2f_kick%.3f_alpha%.2f.txt' % (
                L, steps, T, n_iter, n_periods, J[i], V, hz, kick, alpha[j])

            print(filename[q, i, j])
            files[q, i, j] = np.genfromtxt(filename[q, i, j], skip_header=9, skip_footer=2)
            print(files[q, i, j])
            data[q, i, j] = files[q, i, j]
            print(data[q, i, j][0, :])

#"""

for i in range(num_J):
    for j in range(num_alpha):
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])
            tau[q, i, j] = data[q, i, j][:, 1]


tau_avg = np.zeros((numL, num_J, num_alpha))
tau_err = np.zeros((numL, num_J, num_alpha))


for q in range(numL):
    L = Larray[q]
    for i in range(num_J):
        for j in range(num_alpha):
            tau_avg[q, i, j] = np.mean(tau[q, i, j])
            tau_err[q, i, j] = np.std(tau[q, i, j]) / sqrt(np.size(tau[q,i,j]))


cmJ = cm.get_cmap('turbo', num_J)
    
# Plot tau vs L with error bars
for j in range(num_alpha):
    m += 1
    plt.figure(m)
    for i in range(num_J):
        print([J[i], alpha[j]])
        
        string = '$J = %.g$' % J[i]
        plt.errorbar(Larray[:] , tau_avg[:, i, j], tau_err[:, i, j],
                     label=string, linewidth=1, color=to_hex(cmJ(i)), marker='o')
                     
    plt.xticks(Larray)
    plt.yscale('log')
    plt.xlabel('$L$', fontsize=12)
    plt.ylabel('$\\overline{\\tau}$', fontsize=12)
    plt.legend(fontsize=9, loc='right', bbox_to_anchor=(1.3,0.5))
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)
    plt.ylim([10**1, 10**9])
    
    txt = 'figures/Decay_Times_avg_wrt_L_kick%.2f_alpha%.2f.pdf' % (kick, alpha[j])
    print(txt)
    plt.savefig(txt, dpi=600, bbox_inches='tight')



X = np.empty((numL, num_J, num_alpha), dtype=object)


X = tau

for i in range(num_J):
    for j in range(num_alpha):
        m += 1
        plt.figure(m)
        minbin = min( np.array([ min(log10(X[q,i,j])) for q in range(numL) ], dtype=float) )
        maxbin = max( np.array([ max(log10(X[q,i,j])) for q in range(numL) ], dtype=float) )
        binwidth = 0.5
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])

            string = '$L = %d$' % L
            weights = np.ones_like(X[q, i, j])/len(X[q, i, j])
            plt.hist(X[q, i, j], weights=weights, bins=10**np.arange(minbin, maxbin + binwidth, binwidth),
                     linewidth=1.2, label=string, alpha=0.5)
            #plt.hist(X[q, i, j], weights=weights, linewidth=1.2, label=string, alpha=0.5)

        plt.xlabel('$\\overline{\\tau}$', fontsize=12)
        plt.ylabel('$P(\\overline{\\tau})$', fontsize=12)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.xscale('log')
        #plt.ylim([0, 1])
        #plt.xlim([-18,1])
        plt.legend()

        txt = 'figures/Hist_Decay_Times_wrt_L_J%.5f_kick%.2f_alpha%.2f.pdf' % (J[i], kick, alpha[j])
        print(txt)
        plt.savefig(txt, dpi=600, bbox_inches='tight')
        #plt.close()
        print("\n")
#"""
