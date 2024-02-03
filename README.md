# Ising_Model_Julia

Ising model, a simple statistical physics model that exhibits phase transition, was first proposed by William Lenz in 1920, and its one dimensional case was exactly solved
in 1924 by his student Ernst Ising (however, there is no phase transition in 1D Ising model) [1]. Phase transition for higher dimensional Ising models has been extensively
studied, and exact analytical solution for 2D Ising model in the absence of external magnetic field was found by Lars Onsager in 1944 [1].

An example of phase transition where Ising model has been widely used is magnetism. In a magnetic material, each of the atoms have an associated magnetic spin. Depending on whether the spins are aligned or randomly distributed, the material can have a net magnetic moment (ferromagnetic) or no magnetic moment (paramagnetic). In the absence of an external magnetic field, paramagnetism is observed at high temperatures and ferromagnetism at low temperature. The critical temperature for this magnetic phase transition is called Curie temperature (T_c).

In this project, phase transition for 2D Ising model was studied numerically using Monte Carlo simulations performed via two different algorithms, namely Metropolis and Wolff. Autocorrelation times, critical slowing down phenomena, spin-spin correlation function and finite size scaling were investigated for the 2D lattice system. Phase transition characteristics for a 3D Ising model were also evaluated. All the simulations and analysis was performed in Julia.

References:
[1] Harvey Gould and Jan Tobochnik. Magnetic systems. In Statistical and Thermal Physics, pages 231–293. Princeton University Press, 2010.
