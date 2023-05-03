# Swap

This project deals with computing the evolution and spectrum of periodically driven (Floquet) spin chains; also, as a by product, with Many-Body Localization (MBL) spin models, for both the spin-1/2 and spin-1 cases.
Here we generate Spin Hamiltonians for disordered systems and short-range interactions, both as a Dense Matrix and a Sparse Matrix and study the resulting dynamics and spectral properties. We usually consider hamiltonians with a Z magnetic field and XXZ Heisenberg interactions:

$$ H = \sum_k h_k\sigma_k^z + V_k\sigma_k^z\sigma_{k+1}^z + J_k(\sigma_k^x\sigma_{k+1}^x + \sigma_k^y\sigma_{k+1}^y). $$

Using Exact Diagonalization (ED), we can compute the entire spectrum, which allows us to compute the time evolution of an initial (pure) state by computing the matrix exponential of the Hamiltonian, and apply it iteratively to find the dynamics.

$$ H = WEW^\dagger \Longrightarrow e^{-iH} = We^{-iE}W^\dagger, $$

$$ |\psi((n+1)T)\rangle = U_F|\psi(nT)\rangle. $$

We can also study the spectral properties (distribution of level spacing, average gap ratio, entanglement in eigenstates and more), which have been implemented from scratch, using basic Fortran.

Using Krylov-based methods (through the EXPOKIT library), we can also compute the evolution for larger spin chains, since the Hamiltonians are very sparse (short-range off-diagonal interactions). The resulting dynamics takes however much longer to compute.

Finally, when there are conserved quantities (e.g. \$S_z\$), we can restrict the analysis to a given subspace, so as to further improve performance.

The main (dynamic) observable studied in the evolution is the local magnetization $\langle \sigma_k^z(t)\rangle $, for which we compute the disorder- and time-average, so as to obtain finite-size scaling plots which point to a (continuous) phase transition. We usually deal wih global quantities, such as the imbalance:

$$I = \sum_k (-1)^k \sigma_k^z.$$

