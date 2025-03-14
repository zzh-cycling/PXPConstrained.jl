using PXPConstrained
using LinearAlgebra
using BenchmarkTools
norm(rotated_psi_state_mss(12,0,0.0))
s = rotated_psi_state_mss(12,0,π/3)
norm(s)





norm(rotated_psi_state(12,0.0))
rotated_psi_state(12,π)
using BenchmarkTools
norm(rotated_psi_state(12,π))
a = findfirst(x->x≈1,rotated_psi_state(12,π))
print
norm(rotated_psi_state_mss(12,0,0.0))