include("Ising_2D_wolff.jl")
using .Ising_2D_wolff

N_sweeps = 15000
N_skip = 100
N_thermalize = 3000
L = 100
init = "uniform"
T = 2.3
corr_check = false

m, x, cv = Ising_2D_wolff.main(N_sweeps, N_skip, N_thermalize, L, init, T, corr_check)