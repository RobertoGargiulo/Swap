import numpy as np
import matplotlib.pyplot as plt
from scipy.special import comb
from matplotlib.colors import to_hex
import matplotlib.cm as cm
#from matplotlib.backends.backend_pdf import PdfPages
from numpy import sqrt

print("-----------------------------------")

# Set formatting options
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams['figure.figsize'] = [5,3]

# Import Files
Larray = np.arange(4, 11, 2)
T = 1
#J = 0.005 * 2 ** np.arange(0, 11)
J = np.array([0.01, 0.04, 0.16, 0.32, 0.64, 1.28 ])
V = 3
hz = 16
kick = 0.00
alpha = [0.50, 3.00]
num_J = len(J)
num_alpha = len(alpha)
n_points = num_J * num_alpha
numL = len(Larray)
iter2 = 5120

m = 0

filename = np.empty((numL, num_J, num_alpha), dtype=object)
files = np.empty((numL, num_J, num_alpha), dtype=object)
data = np.empty((numL, num_J, num_alpha), dtype=object)
arr = np.empty((numL, 1), dtype=object)

CORR = np.empty((numL, num_J, num_alpha), dtype=object)
CORR2 = np.empty((numL, num_J, num_alpha), dtype=object)
LI = np.empty((numL, num_J, num_alpha), dtype=object)
IPR = np.empty((numL, num_J, num_alpha), dtype=object)
QE = np.empty((numL, num_J, num_alpha), dtype=object)


edges = np.empty((numL, num_J, num_alpha), dtype=object)
counts = np.empty((numL, num_J, num_alpha), dtype=object)

for i in range(num_J):
    for j in range(num_alpha):
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])
            n_iter = iter2 * 2 ** (1 - L / 2)

            filename[q, i, j] = 'data/eigen/correlations_Swap_LR_Sz0_nspin%d_period%.2f_n_disorder%d_Jxy%.5f_Vzz%.2f_hz%.2f_kick%.3f_alpha%.2f.txt' % (
                L, T, n_iter, J[i], V, hz, kick, alpha[j])

            print(filename[q, i, j])
            files[q, i, j] = np.genfromtxt(filename[q, i, j], skip_header=8)
            print(files[q, i, j])
            data[q, i, j] = files[q, i, j]
            print(data[q, i, j][0, :])

for i in range(num_J):
    for j in range(num_alpha):
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])
            n_iter = iter2 * 2 ** (1 - L / 2)
            CORR[q, i, j] = data[q, i, j][:, 1]
            CORR2[q, i, j] = data[q, i, j][:, 2]
            LI[q, i, j] = data[q, i, j][:, 3]
            IPR[q, i, j] = data[q, i, j][:, 4]
            #QE[q, i, j] = data[q, i, j][:, 5]


# Calculate averages and errors
IPRavg = np.zeros((numL, num_J, num_alpha))
IPRerr = np.zeros((numL, num_J, num_alpha))
CORRavg = np.zeros((numL, num_J, num_alpha))
CORRerr = np.zeros((numL, num_J, num_alpha))
CORR2avg = np.zeros((numL, num_J, num_alpha))
CORR2err = np.zeros((numL, num_J, num_alpha))

n_iter = np.array([ iter2 * 2 ** (1 - L / 2) for L in Larray])

for q in range(numL):
    L = Larray[q]
    n_iter = int(iter2 * 2 ** (1 - L / 2))
    for i in range(num_J):
        for j in range(num_alpha):
            IPRavg[q, i, j] = -np.mean(np.log(IPR[q, i, j]))
            IPRerr[q, i, j] = np.std(np.log(IPR[q, i, j]))
            norm = L
            CORRavg[q, i, j] = np.mean(CORR[q, i, j] / norm)
            CORRerr[q, i, j] = np.std(CORR[q, i, j] / norm) / sqrt(np.size(CORR[q,i,j]))
            CORR2avg[q, i, j] = np.mean(CORR2[q, i, j] / norm)
            CORR2err[q, i, j] = np.std(CORR2[q, i, j] / norm) / sqrt(np.size(CORR2[q,i,j]))

cmJ = cm.get_cmap('turbo', num_J)

# Plot CORR vs L with error bars
for j in range(num_alpha):
    m += 1
    plt.figure(m)
    for i in range(num_J):
        print([J[i], alpha[j]])

        string = '$J = %.2f$' % J[i]
        plt.errorbar(Larray[:] , CORRavg[:, i, j], CORRerr[:, i, j],
                     label=string, linewidth=1, color=to_hex(cmJ(i)), marker='o')

    plt.xticks(Larray)
    plt.xlabel('$L$', fontsize=12)
    plt.ylabel('$\overline{\Sigma}$', fontsize=12)
    plt.legend(fontsize=10, loc='lower center', bbox_to_anchor=(0.5,-0.4), fancybox=True, shadow=True, ncol=3)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    txt = 'figures/CORR_sigmaz_avg_wrt_J_kick%.2f_alpha%.2f.pdf' % (kick, alpha[j])
    print(txt)
    plt.savefig(txt, dpi=600, bbox_inches='tight')

for j in range(num_alpha):
    m += 1
    plt.figure(m)
    for i in range(num_J):
        print([J[i], alpha[j]])

        string = '$J = %.3f$' % J[i]
        plt.errorbar(Larray[:] , CORR2avg[:, i, j], CORR2err[:, i, j],
                     label=string, linewidth=1, color=to_hex(cmJ(i)), marker='o')

    plt.xticks(Larray)
    plt.xlabel('$L$', fontsize=12)
    plt.ylabel('$\overline{\Sigma}_2$', fontsize=12)
    plt.legend(fontsize=10, loc='lower center', bbox_to_anchor=(0.5,-0.4), fancybox=True, shadow=True, ncol=3)
    plt.xticks(fontsize=10)
    plt.yticks(fontsize=10)

    txt = 'figures/CORR_sigmaz2_avg_wrt_J_kick%.2f_alpha%.2f.pdf' % (kick, alpha[j])
    print(txt)
    plt.savefig(txt, dpi=600, bbox_inches='tight')

X = np.empty((numL, num_J, num_alpha), dtype=object)



for i in range(num_J):
    for j in range(num_alpha):
        m += 1
        plt.figure(m)
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])

            X[q, i, j] = CORR[q, i, j]
            string = '$L = %d$' % L
            weights = np.ones_like(X[q, i, j])/len(X[q, i, j])
            binwidth = 1.5
            plt.hist(X[q, i, j], weights=weights, bins=np.arange(min(X[q,i,j]), max(X[q,i,j]) + binwidth, binwidth),
                     linewidth=1.2, label=string, alpha=0.5)

        plt.xlabel('$\\Sigma_{\\alpha}$', fontsize=12)
        plt.ylabel('$P(\\Sigma_{\\alpha})$', fontsize=12)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend()

        txt = 'figures/Hist_Correlations_sigmaz_wrt_L_J%.5f_kick%.2f_alpha%.2f.pdf' % (J[i], kick, alpha[j])
        print(txt)
        plt.savefig(txt, dpi=600, bbox_inches='tight')
        print("\n")

for i in range(num_J):
    for j in range(num_alpha):
        m += 1
        plt.figure(m)
        for q in range(numL):
            L = Larray[q]
            print([J[i], alpha[j], L])

            X[q, i, j] = CORR2[q, i, j]
            string = '$L = %d$' % L
            weights = np.ones_like(X[q, i, j])/len(X[q, i, j])
            binwidth = 1.5
            plt.hist(X[q, i, j], weights=weights, bins=np.arange(min(X[q,i,j]), max(X[q,i,j]) + binwidth, binwidth),
                     linewidth=1.2, label=string, alpha=0.5)

        plt.xlabel('$\\Sigma_{2,\\alpha}$', fontsize=12)
        plt.ylabel('$P(\\Sigma_{2,\\alpha})$', fontsize=12)
        plt.xticks(fontsize=10)
        plt.yticks(fontsize=10)
        plt.legend()

        txt = 'Figures/SpinHalf/Hist_Correlations_sigmaz2_wrt_L_J%.5f_kick%.2f_alpha%.2f.pdf' % (J[i], kick, alpha[j])
        print(txt)
        plt.savefig(txt, dpi=600, bbox_inches='tight')
        print("\n")
