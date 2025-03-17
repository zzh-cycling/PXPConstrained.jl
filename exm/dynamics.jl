using PXPConstrained
using LinearAlgebra
using BenchmarkTools
a=norm(rotated_psi_state_mss(12,0,0.0))
s = rotated_psi_state_mss(12,0,π/3)
norm(s)
q = rotated_psi_state_mss(12,0,π)
norm(q)



norm(rotated_psi_state(12,0.0))
rotated_psi_state(12,π)
norm(rotated_psi_state(12,π))
norm(rotated_psi_state(12,π/2))
norm(rotated_psi_state_mss(12,0,0.0))


using BitBasis
function rotated_psi_state_mss2(::Type{T}, k::Int64, θ::Real) where {N, T<: BitStr{N}}
    # params: a state in maximum symmetry space, and the momentum of the state
    # return: the state in total space
    MSS, MSS_dic, qlist = PXP_MSS_basis(T, k)
    γ = tan(θ/2)
    basisK, k_dic = PXP_K_basis(T, k)
    
    rotated_state = zeros(Float64, length(MSS))

    for (i, base) in enumerate(MSS)
        Y= sqrt(length(k_dic[base]))/N
        Z= sqrt(qlist[i])*Y/2
        amp1 = even_zeros(base, γ)
        amp2 = even_ones(base, γ)
        
        rotated_state[i]=Z*N*sqrt(2)*(amp1+amp2)
    end

    # rotated_state = rotated_state .* cos(θ/2)^N
    
    return rotated_state /=  norm(rotated_state .* cos(θ/2)^N)
end
rotated_psi_state_mss2(N::Int64, k::Int64, θ::Real) = rotated_psi_state_mss2(BitStr{N, Int}, k, θ)
o = rotated_psi_state_mss2(12,0,π/2)
norm(o)





