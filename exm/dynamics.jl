using PXPConstrained
using LinearAlgebra
using BenchmarkTools
a=norm(rotated_psi_state_mss(12,0,0.0))
s = rotated_psi_state_mss(12,0,π/3)
norm(s)





norm(rotated_psi_state(12,0.0))
rotated_psi_state(12,π)
norm(rotated_psi_state(12,π))
norm(rotated_psi_state_mss(12,0,0.0))