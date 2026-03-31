using PXPConstrained
using LinearAlgebra
using JLD2, JLD

lambdalis1=collect(0:0.1:1)
lambdalis2=collect(0.25:0.01:0.75)
lambdalis=vcat(lambdalis1[1:3],lambdalis2,lambdalis1[end-2:end])

function load_data_m(N)
    data = load("/hpc2hdd/home/zzhi359/big_data/PXP_Ham_wf/eigen_MSS_N$(N).jld2")
    energy, states = data["energy"], data["states"]

    index_table = Dict(
        18 => 136,
        22 => 655,
        26 => 3502,
        30 => 19934,
    )
    dE_table = Dict(
        18 => 5e-1,
        22 => 1e-1,
        26 => 2e-2,
        30 => 4e-3,
    )
    scar = states[:, index_table[N]]
    dE = dE_table[N]
    thermal_inwindow_index = setdiff(findall(x -> abs(x - energy[index_table[N]]) < dE, energy), [index_table[N]])
    thermal_states = states[:, thermal_inwindow_index]

    @show [scar'*i for i in eachcol(thermal_states)]
    return scar, thermal_states
end

function load_data_p(N)
    data = load("/hpc2hdd/home/zzhi359/big_data/PXP_Ham_wf/eigen_MSS_N$(N).jld2")
    energy, states = data["energy"], data["states"]

    index_table = Dict(
        22 => 655,
        26 => 3502,
        30 => 9913,
    )
    dE_table = Dict(
        22 => 1e-1,
        26 => 2e-2,
        30 => 1e-3,
    )
    scar = states[:, index_table[N]]
    dE = dE_table[N]
    thermal_inwindow_index = setdiff(findall(x -> abs(x - energy[index_table[N]]) < dE, energy), [index_table[N]])
    thermal_states = states[:, thermal_inwindow_index]

    # Inv_op = inversion_matrix(N)
    # Trans_op = translation_matrix(N)
    # l = size(Trans_op, 2)

    # mom_proj = momentum_projector(Trans_op, N; kpi=false)
    # inv_proj = (Matrix{Float64}(I, l, l) + Inv_op) / 2
    # sector_proj = inv_proj * mom_proj

    # scar = sector_proj * scar
    # scar /= norm(scar)

    # # 注意：投影不是幺正变换，正交态投影后一般不会保持两两正交，需要再做一次正交化。
    # thermal_states_symmetry = project_and_orthonormalize(thermal_states, sector_proj)
    @show [scar'*i for i in eachcol(thermal_states)]
    return scar, thermal_states
end

function superposition(scar::Vector{Float64}, thermal_ensemble::Matrix{Float64}, lambda::Float64)
    state = (1 .- lambda) .* scar .+ lambda .* thermal_ensemble
    norms =  mapslices(norm, state; dims=1)
    state=state./norms
    # return all the superposition states with different lambda, which is a matrix with each column being a state.
    return state
end # equivalent to superposite seperately.

function WA_den_scaling(N::Int64, l::Int64, super_states::Matrix{Float64}) # NOT Density
    HA = PXP_Ham(l, false)
    subenergy = eigvals(HA)[1]
    GS_energy=0.0
    passive_energy=0.0
    len = size(super_states)[2]

    for j in 1:len
        @show j, "WA"
        state = super_states[:, j]
        gs, sub, pass = ergotropy_PXP_MSS_state(N, l, state, 0)
        GS_energy+=gs
        passive_energy+=pass
        # linear average over the superposition states
    end

    GS_energy = GS_energy / len
    passive_energy = passive_energy / len

    WA_den = (GS_energy-passive_energy)/N
    QA_den = (passive_energy - subenergy[1])/N
    DeltaEA_den = (GS_energy - subenergy[1])/N
    return WA_den, QA_den, DeltaEA_den
end

function EE_density_scaling(N::Int64, super_states::Matrix{Float64})
    EE = 0.0
    # half chain entanglement entropy density
    len = size(super_states)[2]

    for j in 1:len
        @show j, "EE"
        state = super_states[:, j]
        subrho = rdm_PXP_MSS(N, collect(1:div(N,2)), state, 0)
        EE+=ee(subrho)
    end

    EE = EE / len
    EE_den = EE/N
    return EE_den
end

function ergo_scaling(N)
    data_path = joinpath("data/transition_symmetry_super/new_scar") 
    save_path = joinpath("data/transition_symmetry_super/new_scar")

    WA_lis=zeros(length(lambdalis))
    QA_lis=zeros(length(lambdalis))
    DeltaEA_lis=zeros(length(lambdalis))
    EE_den_lis=zeros(length(lambdalis))

    # Because of the symmetry sector we choose, the next scar state is in the minus inversion symmetry sector for N=12,16,20,22, so we need to load the data with minus inversion symmetry.
    if N % 4 == 0
        scar, thermal_ensemble = load_data_p(N)
    else
        scar, thermal_ensemble = load_data_m(N)
    end
    
    for index in eachindex(lambdalis)
        super_states = superposition(scar, thermal_ensemble, lambdalis[index])
        WA, QA, DeltaEA = WA_den_scaling(N, div(N,2), super_states)
        WA_lis[index]=WA
        QA_lis[index]=QA
        DeltaEA_lis[index]=DeltaEA
        EE_den_lis[index] = EE_density_scaling(N, super_states)
    end
    save(joinpath(save_path, "ergo_denlis_N$(N).jld"), "ergo_den_lis", WA_lis, "boundenergy_den_lis", QA_lis, "excessenergy_den_lis", DeltaEA_lis, "ee_den_lis", EE_den_lis)
end

if length(ARGS) == 0
    println("No arguments provided.")
else
    arg=parse.(Int64, ARGS)[1]
    println("Received argument: N=$arg")
    ergo_scaling(arg)
end