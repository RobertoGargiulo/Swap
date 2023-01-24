# Swap

This project deals with computing the evolution and spectrum of periodically driven (Floquet) spin chains; also, as a by product, with Many-Body Localization (MBL) spin models, for both the spin-1/2 and spin-1 cases.
Here we generate Spin Hamiltonians for disordered systems and short-range interactions, both as a Dense Matrix and a Sparse Matrix and study the resulting dynamics and spectral properties.\

Using Exact Diagonalization (ED), we can compute the entire spectrum, which allows us to compute the time evolution of an initial (pure) state by computing the matrix exponential of the Hamiltonian, and apply it iteratively to find the dynamics. We can also study the spectral properties (distribution of level spacing, average gap ratio, entanglement in eigenstates and more).

Using Krylov-based methods, we can also compute the evolution for larger spin chains, since the Hamiltonians are very sparse (short-range interactions). The resulting dynamics takes however much longer to compute.

Finally, when there are conserved quantities (e.g. \$S_z\$), we can restrict the analysis to a given subspace, so as to further improve performance.

The main (dynamic) observable studied in the evolution is the imbalance, for which we compute the disorder- and time-average, so as to obtain finite-size scaling plots which point to a (continuous) phase transition.

