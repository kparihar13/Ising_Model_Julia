include("Ising_2D_metropolis.jl")
using .Ising_2D_metropolis

N_sweeps = 5000
N_skip = 100
N_thermalize = 500
L = 100
init = "uniform"
T = 2.0
corr_check = false

m, x, cv = Ising_2D_metropolis.main(N_sweeps, N_skip, N_thermalize, L, init, T, corr_check)