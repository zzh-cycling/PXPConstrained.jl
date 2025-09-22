"""
    Dynamics.jl

Functions for quantum time evolution and special local rotated state preparation in the PXP model.
This module provides tools for computing time evolution using both exact diagonalization
and sparse matrix methods, as well as utilities for preparing special quantum states.
"""

"""
    wf_time_evolution(psi0::Vector{T}, times::Vector{Float64}, energy::Vector{Float64}, states::Matrix{Float64}) where {T <: Real}

Compute time evolution of a quantum state using eigendecomposition.

Evolves the initial state |ψ(0)⟩ according to |ψ(t)⟩ = e^{-iHt}|ψ(0)⟩
using the eigenbasis decomposition. Uses multi-threading for parallel computation.

# Arguments
- `psi0::Vector{T}`: Initial state vector
- `times::Vector{Float64}`: Time points to evaluate
- `energy::Vector{Float64}`: Eigenvalues of the Hamiltonian
- `states::Matrix{Float64}`: Eigenvectors of the Hamiltonian (columns)

# Returns
- `Vector{Vector{ComplexF64}}`: Time-evolved state at each time point

# Example
```julia
H = PXP_Ham(10)
energy, states = eigen(H)
psi0 = zeros(length(PXP_basis(10)))  # all zero state
times = collect(0:0.1:10)
psi_t = wf_time_evolution(psi0, times, energy, states)
```
"""
function wf_time_evolution(psi0::Vector{T}, times::Vector{Float64}, energy::Vector{Float64},states::Matrix{Float64}) where {T <: Real}
    wflis=Vector{Vector{ComplexF64}}(undef,length(times))
    c = states'*psi0
    exp_factors = [exp.(-1im * t * energy) for t in times]
    
    # Use multi-threading for parallel computation
    Threads.@threads for i in eachindex(times)
        wflis[i] = states * (c .* exp_factors[i])
    end
    return wflis
end

"""
    wf_time_evolution_sparse(L::Int, k::Int64, psi0::Vector{ET}, t_values::Vector{Float64}) where {ET}

Compute time evolution using sparse matrix exponentiation.

Efficiently evolves quantum states in the maximum symmetry subspace using
the Krylov subspace method for matrix exponential action. More memory-efficient
than full diagonalization for large systems.

# Arguments
- `L::Int`: System size
- `k::Int64`: Momentum quantum number
- `psi0::Vector{ET}`: Initial state in MSS basis
- `t_values::Vector{Float64}`: Time points (assumed equally spaced)

# Returns
- `Vector{Vector{ComplexF64}}`: Time-evolved states at each time point

# Example
```julia
# Initial Z2 state in MSS basis
using PXPConstrained, BitBasis
psi0_mss = zeros(length(PXP_MSS_basis(BitStr{12, Int}, 0)[1]))
psi0_mss[end] = 1  # Example state in MSS basis
times = collect(0:0.01:5.0)
psi_t = wf_time_evolution_sparse(12, 0, psi0_mss, times)
```
"""
function wf_time_evolution_sparse(L::Int, k::Int64, psi0::Vector{ET}, t_values::Vector{Float64}) where {ET}
    # 预先计算时间步长
    dt = t_values[2] - t_values[1]
    H = PXP_MSS_Ham_sparse(L, k)
    # 初始化存储波函数的数组
    @assert length(psi0) == size(H, 1) "state length is expected to be $(size(H, 1)), but got $(length(psi0))"
    wavefunctions = Vector{Vector{ComplexF64}}(undef, length(t_values))
    wavefunctions[1] = copy(psi0)
        
    initial_norm = norm(psi0)
    if !isapprox(initial_norm, 1.0, atol=1e-10)
        @warn "The initial state isn't normalized, with norm: $initial_norm"
        psi0 = psi0 / initial_norm
    end
        
    psi_t = copy(psi0)
    for i in 2:length(t_values)
        # 使用 expv 计算波函数的时间演化
        psi_t = expv(-im * dt, H, psi_t; ishermitian=true)   
        wavefunctions[i] = psi_t
    end
    
    return wavefunctions
end

"""
    rotated_psi_state(::Type{T}, θ::Real) where {N, T<: BitStr{N}}
    rotated_psi_state(N::Int64, θ::Real)

Generate a rotated Z2 quantum state in the PXP basis.

Creates the state obtained by applying the rotation exp(iθ/2 Y) to the
Z2 state |10101010...⟩, where Y is the total spin-Y operator.

# Arguments
- `T::Type{BitStr{N}}` or `N::Int64`: System size specification  
- `θ::Real`: Rotation angle

# Returns
- `Vector{Float64}`: Normalized rotated state in PXP basis

# Example
```julia
# Create rotated Z2 state at π/3
psi_rot = rotated_psi_state(10, π/3)
```
"""
function rotated_psi_state(::Type{T}, θ::Real, pbc::Bool=true) where {N, T<: BitStr{N}}
    # params: the particlenumber of the space, and rotation angle θ for the Z2 state
    # return: the state rotated by on site rotation exp(i θ/2 Y)
    basis = PXP_basis(T, pbc)
    rotated_state = zeros(Float64, length(basis))
    
    for (i, base) in enumerate(basis)
        exp = count_zeros_and_ones(base)
        rotated_state[i] = Z2_overlap(exp, θ)
    end
    
    return rotated_state ./ norm(rotated_state)
end
rotated_psi_state(N::Int64, θ::Real, pbc::Bool=true) = rotated_psi_state(BitStr{N, Int}, θ, pbc)


function count_zeros_and_ones(base::BitStr{N}) where {N}
    even_zeros = 0
    even_ones = 0
    odd_zeros = 0
    odd_ones = 0
    for j in 1:N
        if j % 2 == 1  # 奇数位
            if base[j] == 1
                odd_ones += 1
            else
                odd_zeros += 1
            end
        else  # 偶数位
            if base[j] == 1
                even_ones += 1
            else
                even_zeros += 1
            end
        end
    end

    return even_zeros, even_ones, odd_zeros, odd_ones
end


"""
    Z2_overlap(exp::Tuple{Int64, Int64, Int64, Int64}, θ::Real)

Compute overlap amplitude for rotated Z2 state |10101010...⟩.

Calculates the amplitude ⟨basis|e^{iθ/2 Y}|10101010...⟩ where Y is the
total spin-Y operator. Used for constructing rotated Z2 states.

# Arguments
- `exp::Tuple{Int64, Int64, Int64, Int64}`: (even_zeros, even_ones, odd_zeros, odd_ones)
- `θ::Real`: Rotation angle

# Returns
- `Float64`: Matrix element amplitude

# Example
```julia
counts = (1, 1, 1, 1)
amplitude = Z2_overlap(counts, π/4)
```
"""
#对于｜10101010...>态，做rotation后的振幅计算,从右往左计数
function Z2_overlap(exp::Tuple{Int64, Int64, Int64, Int64}, θ::Real)
    even_zeros, even_ones, odd_zeros, odd_ones = exp
    return sin(θ/2)^(odd_ones) * (-sin(θ/2))^even_zeros * cos(θ/2)^(even_ones+odd_zeros)
end

"""
    Z2tilde_overlap(exp::Tuple{Int64, Int64, Int64, Int64}, θ::Real)

Compute overlap amplitude for rotated Z2 state |01010101...⟩.

Calculates the amplitude ⟨basis|e^{iθ/2 Y}|01010101...⟩ where Y is the
total spin-Y operator. Used for constructing rotated anti-Z2 states.

# Arguments
- `exp::Tuple{Int64, Int64, Int64, Int64}`: (even_zeros, even_ones, odd_zeros, odd_ones)
- `θ::Real`: Rotation angle

# Returns
- `Float64`: Matrix element amplitude

# Example
```julia
counts = (1, 1, 1, 1)
amplitude = Z2tilde_overlap(counts, π/4)
```
"""
#对于｜01010101...>态，做rotation后的振幅计算，从右往左计数
function Z2tilde_overlap(exp::Tuple{Int64, Int64, Int64, Int64}, θ::Real)
    even_zeros, even_ones, odd_zeros, odd_ones = exp
    return sin(θ/2)^(even_ones) * (-sin(θ/2))^odd_zeros * cos(θ/2)^(even_zeros+odd_ones)
end

function rotated_psi_state_mss(::Type{T}, k::Int64, θ::Real, inv::Int64=1) where {N, T<: BitStr{N}}
    # params: a state in maximum symmetry space, and the momentum of the state
    # return: the state in total space
    @assert k == 0 || k==div(N,2) "k is expected to be 0 or $(div(N,2)), but got $k"
    @assert inv ==1 || inv==-1 "inv is expected to be 1 or -1, but got $(inv)"
    MSS, MSS_dic, qlist = PXP_MSS_basis(T, k, inv)
    basisK, k_dic = PXP_K_basis(T, k)
    
    rotated_state = zeros(Float64, length(MSS))

    if inv==1 && k==0 || inv==-1 && k==div(N,2)
        for (i, base) in enumerate(MSS)
            Y = sqrt(length(k_dic[base]))/N
            Z = sqrt(qlist[i])*Y/2
            
            # 计算基态在旋转后的振幅
            exp = count_zeros_and_ones(base)
            amp1 = Z2_overlap(exp, θ)
            amp2 = Z2tilde_overlap(exp, θ)
            rotated_state[i] = Z*N*sqrt(2)*(amp1+amp2)
            
        end
    elseif inv==-1 && k==0 || inv==1 && k==div(N,2)
        for (i, base) in enumerate(MSS)
            Y = sqrt(length(k_dic[base]))/N
            Z = sqrt(qlist[i])*Y/2
            
            # 计算基态在旋转后的振幅
            exp = count_zeros_and_ones(base)
            amp1 = Z2_overlap(exp, θ)
            amp2 = Z2tilde_overlap(exp, θ)
            rotated_state[i] = Z*N*sqrt(2)*(amp1-amp2)
            
        end
    end

    # return rotated_state
    return rotated_state ./ norm(rotated_state)
end
rotated_psi_state_mss(N::Int64, k::Int64, θ::Real, inv::Int64=1) = rotated_psi_state_mss(BitStr{N, Int}, k, θ, inv)
