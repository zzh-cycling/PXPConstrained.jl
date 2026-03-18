using PXPConstrained
using LinearAlgebra

scar_indexlis16=[1, 2, 9, 27, 82, 202, 408, 728, 1075, 1480, 1800, 2006, 2126, 2181, 2199, 2206, 2207]

H = PXP_Ham(16)
energy, states = eigen(H)

overlaplis  = log.(states[end,:].^2)
scar_overlap = overlaplis[scar_indexlis16]

zero_energy_indices = findall(x->isapprox(x, 0.0, atol=1e-10), energy)

ΔE1 = energy[scar_indexlis16[8]-1] - energy[scar_indexlis16[8]]
ΔE2 = energy[scar_indexlis16[8]-2] - energy[scar_indexlis16[8]]
ΔEp1 = energy[scar_indexlis16[8]+1] - energy[scar_indexlis16[8]]
ΔEp2 = energy[scar_indexlis16[8]+2] - energy[scar_indexlis16[8]]

energy_shell_index = findall(x->isapprox(x, energy[scar_indexlis16[8]], atol=0.004), energy)


energy[zero_energy_indices]
energy[zero_energy_indices[1]-1]
energy[zero_energy_indices[end]+1]


H2 = PXP_Ham(18)
energy2, states2 = eigen(H2)

zero_energy_indices2 = findall(x->isapprox(x, 0.0, atol=1e-10), energy2)

energy2[zero_energy_indices2[1]]
energy2[zero_energy_indices2[1]-1]
energy2[zero_energy_indices2[end]+1]