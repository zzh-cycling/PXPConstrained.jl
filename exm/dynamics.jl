using PXPConstrained
using LinearAlgebra
using BenchmarkTools
a=norm(rotated_psi_state_mss(12,0,0.0))
rotated_psi_state_mss(12,0,π)

s = rotated_psi_state_mss(12,0,π/3)
norm(s)
q = rotated_psi_state_mss(12,0,0.0)
norm(q)
@show q


norm(rotated_psi_state(12,0.0))
rotated_psi_state(12,π)
norm(rotated_psi_state(12,π))
norm(rotated_psi_state(12,π/2))
norm(rotated_psi_state(12,π/7))

rotated_psi_state_mss(12,0,0.0)


norm(rotated_psi_state_mss(12,0,0.0))

rotated_psi_state(12,0.0)



