ITensors.set_warn_order(60)

function proj_FSA(::Type{T}) where {N, T <: BitStr{N}}
#Generate the projection matrix for the FSA subspace, which should equals to identity matrix in FSA subspace.
#Meanwhile not equal to identity matrix in the total Hilbert space. Input N is the size of the system, return the projection matrix.
    energy,states=eigen(PXP_FSA_Ham(T))

    return states*states'
end

function proj_FSA2total(::Type{T}) where {N, T <: BitStr{N}}
    # iso=iso_total2FSA(N)

    parent_path=homedir()
    file_path=joinpath(parent_path, "/Documents/projects/big_data/scar_thermal_FSA/iso_FSA/iso_total2FSA$(N).jld")

    if isfile(file_path)
        iso = load(file_path, "iso")
    else
        iso = iso_total2FSA(T)
    end

    energy,states=eigen(PXP_FSA_Ham(T))
    Proj=iso*states*states'*iso'
    return Proj
end

function sep_scar_FSA(::Type{T}, energy::Vector{Float64},states::Matrix{Float64}) where {N, T <: BitStr{N}}
    indices=[index for (index,value) in enumerate(energy) if abs(value)<=1e-8]
    P_FSA=proj_FSA2total(T)

    iso_0modes=states[:,indices]
    # Here need to note that for different basic MKL/openblas version, the degenrated subspace will generate different states. So we need to use the same kind of Package.
    # P_0modes=iso_0modes*iso_0modes'

    PPP_symmetrized = (iso_0modes'*P_FSA*iso_0modes + (iso_0modes'*P_FSA*iso_0modes)') / 2
    
    vals, vecs = eigen(PPP_symmetrized)
    vecs=iso_0modes*vecs

    scar,thermal=vecs[:,end],vecs[:,1:end-1]
    # scar is already orthogonal to thermal states, so we do not need to do the gram_schmidt process. (Only one scar, as FSA expected)
    return scar,thermal
end
sep_scar_FSA(n::Int64, energy::Vector{Float64}, states::Matrix{Float64}) = sep_scar_FSA(BitStr{n, Int}, energy, states)

function sep_scar_FSA_inv(::Type{T},energy::Vector{Float64},states::Matrix{Float64}) where {N, T <: BitStr{N}}
    indices=[index for (index,value) in enumerate(energy) if abs(value)<=1e-8]
    P_FSA=proj_FSA2total(T)

    iso_0modes=states[:,indices]
    # We note that the inversion symmetry is exceptionally influence the thermal 

    PPP_symmetrized = (iso_0modes'*P_FSA*iso_0modes + (iso_0modes'*P_FSA*iso_0modes)') / 2

    vals, vecs = eigen(PPP_symmetrized)
    vecs=iso_0modes*vecs

    Inv=inversion_matrix(T)
    states=vecs[:,1:end-1]
    total_states=Vector{Float64}[]
    l=size(states)[1]
    
    for i in 1:size(states)[2]
        st=states[:,i]
        if isapprox(st'*Inv*st, 1, atol=1e-6) || isapprox(st'*Inv*st, -1, atol=1e-6)
            push!(total_states, st)
        else
        stp=(I(l)+Inv)/2*st
        stp/=norm(stp)
        stn=(I(l)-Inv)/2*st
        stn/=norm(stn)
        push!(total_states, stp)
        push!(total_states, stn)
        end
        
    end

    total_states = hcat(total_states...)
    return vecs[:,end], total_states
    # scar,thermal=vecs[:,end],vecs[:,1:end-1]

    # return scar,thermal
end
sep_scar_FSA_inv(n::Int64, energy::Vector{Float64}, states::Matrix{Float64}) = sep_scar_FSA_inv(BitStr{n, Int}, energy, states)

function proj_Ob(energy::Vector{Float64},states::Matrix{Float64},Ob::Matrix{Float64})
# Create the operator in the new basis of zero energy subspaces
# Actually we find that the non correct overlap between the eigenstates and Z2states is due to the superposition of the degenrate scar states 
#and thermal states in zero modes subspace. Not mainly due to the inversion symmetry.
    indices=[index for (index,value) in enumerate(energy) if abs(value)<=1e-8]
    iso_0modes=states[:,indices]
    projected_matrix = iso_0modes'*Ob*iso_0modes
    return projected_matrix
end

function sep_scar_Ob(energy::Vector{Float64},states::Matrix{Float64},Ob::Matrix{Float64})

    project_matrix = proj_Ob(energy, states, Ob);
    proj_val, proj_vec = eigen(project_matrix);
    scar = states[:,indices] * proj_vec[:,end];
    thermal_ensemble=states[:,indices] * proj_vec[:,1:end-1];

    return scar, thermal_ensemble
end

function proj_Z2(energy::Vector{Float64}, states::Matrix{Float64})
    # Create the |Z2state><Z2state| in the new basis
    # Actually we find that the non correct overlap between the eigenstates and Z2states is due to the superposition of the degenrate scar states and thermal states in zero modes subspace. Not mainly due to the inversion symmetry.
    indices=[index for (index,value) in enumerate(energy) if abs(value)<=1e-8]
    iso_0modes=states[:,indices]
    projected_matrix = iso_0modes'*Z2*iso_0modes
    projected_matrix = zeros(length(indices), length(indices))
    for i in eachindex(indices)
        for j in i: length(indices)
            u = real(states[:, indices[i]])
            v = real(states[:, indices[j]])
            # Because we are in the constrained subspace, the overlap between the states and Z2state is just the last element of the state. p[i,j]=u[end]*v[end] basically saying 
            #<state i|Z2state><Z2state|state j>
            projected_matrix[i, j] = u[end] * v[end] #if we choose u[1597] * v[1597] we get the same result, not only overlap and entanglement entropy.
            projected_matrix[j, i] = projected_matrix[i, j]
        end
    end
    return projected_matrix
end


function proj_invZ2(::Type{T}, states::Matrix{Float64},indices::Vector{Int64}) where {N, T <: BitStr{N}}
# Create the operator in the new basis
# Actually we find that the non correct overlap between the eigenstates and Z2states is due to the superposition of the degenrate scar states and thermal states in zero modes subspace. Not mainly due to the inversion symmetry.
    basis = PXP_basis(T)
    projected_matrix = zeros(length(indices), length(indices))
    for i in eachindex(indices)
        for j in i: length(indices)
            u = real(states[:, indices[i]])
            v = real(states[:, indices[j]])
            # Because we are in the constrained subspace, the overlap between the states and Z2state is just the last element of the state. p[i,j]=u[end]*v[end] basically saying 
            #<state i|Z2state><Z2state|state j>
            z2tilde=findfirst(x->x==div(2^N-1, 3),basis)+1
            projected_matrix[i, j] = u[end] * v[end]+u[z2tilde] * v[z2tilde]+u[end] *v[z2tilde]+ u[z2tilde] *v[end]#if we choose u[1597] * v[1597] we get the same result, not only overlap and entanglement entropy.
            projected_matrix[j, i] = projected_matrix[i, j]
        end
    end
    return projected_matrix
end


function gene_scar(N::Int64)
    B0=Matrix{Float64}([1.0 0.0 0.0 ; 0.0 1.0 0.0])
    B1=Matrix{Float64}(sqrt(2).*[0.0 0.0 0.0 ; 1.0 0.0 1.0])
    C0=Matrix{Float64}([0.0 -1.0; 1.0 0.0 ;0.0 0.0])
    C1=Matrix{Float64}(sqrt(2).*[1.0 0.0; 0.0 0.0 ; -1.0 0.0])

    tensors = []
    b_idx = [Index(2,"b1")]
    for i in 1:N
        if isodd(i)
    		push!(b_idx, Index(3,"b$(i+1)"))
            temp=zeros((2,3,2))
            temp[:,:,1]= B0
            temp[:,:,2]=B1
            push!(tensors, ITensor(temp, b_idx[end-1], b_idx[end], Index(2, "p$(i)")))
        else
    		push!(b_idx, Index(2,"b$(i+1)"))
            temp=zeros((3,2,2))
            temp[:,:,1]= C0
            temp[:,:,2]=C1
            push!(tensors, ITensor(temp, b_idx[end-1], b_idx[end], Index(2, "p$(i)")))
        end
    end

    scar = reduce(*, tensors)
    scar = scar*delta(b_idx[1], b_idx[end])
    
    # Exactly equal to the translated scar state T*scar1
end

function vec2k0pi(::Type{T}, state::Vector{ET}) where {N, T <: BitStr{N}, ET}
    basis = PXP_basis(T)
    l=length(basis)
    Trans=translation_matrix(T)
    op1=I(l)
    op2=I(l)


    T_power = I(l)  
    # 计算 Trans 的幂并同时更新 op1 和 op2
    for i in 1:(N-1)
        T_power *= Trans  # 计算 T^i
        op1 += T_power
        op2 += (-1)^i * T_power
    end

    statek0=op1*state
    statek0/=norm(statek0)
    statekpi=op2*state
    statekpi/=norm(statekpi)
    return statek0,statekpi
end

function sep_scar_exact(::Type{T}, energy::Vector{Float64}, states::Matrix{Float64}) where {N, T <: BitStr{N}} 
    indices = [index for (index, value) in enumerate(energy) if abs(value) <= 1e-8]
    Trans = translation_matrix(T)
    scar=gene_scar(N)
    scar1=storage(scar)

    basis= PXP_basis(T)
    basis_int = [i.buf for i in basis]
    scar1=scar1[basis_int.+1]
    scar1/=norm(scar1)

    scar2=Trans*scar1
    scar2/=norm(scar2)

    exact_scar= scar1+scar2
    exact_scar/=norm(exact_scar)

    exact_scar_prime=scar1-scar2
    exact_scar_prime/=norm(exact_scar_prime)

    Q, R = qr(hcat(exact_scar,exact_scar_prime,states[:,indices]))
    # Noting that not full rank matrix Q needs to take part.
    l=size(R)[1]
    thermal_ensemble=Q[:,3:l]

    return exact_scar, exact_scar_prime, thermal_ensemble
    
end
sep_scar_exact(n::Int64, energy::Vector{Float64}, states::Matrix{Float64}) = sep_scar_exact(BitStr{n, Int}, energy, states)

function sep_scar_exact_translation(::Type{T}, energy::Vector{Float64}, states::Matrix{Float64}) where {N, T <: BitStr{N}}
    indices = [index for (index, value) in enumerate(energy) if abs(value) <= 1e-8]

    Trans = translation_matrix(T)
    scar=gene_scar(N)
    scar1=storage(scar)

    basis= PXP_basis(T)
    basis_int = [i.buf for i in basis]
    scar1=scar1[basis_int.+1]
    scar1/=norm(scar1)

    scar2=T*scar1
    scar2/=norm(scar2)
    
    exact_scar= scar1+scar2
    exact_scar/=norm(exact_scar)

    exact_scar_prime=scar1-scar2
    exact_scar_prime/=norm(exact_scar_prime)

    
    Q, R = qr(hcat(exact_scar,exact_scar_prime,states[:,indices]))
    l=size(R)[1]
    thermal_ensemble=Q[:,3:l]

    total_thek0=Vector{Float64}[]
    total_thekpi=Vector{Float64}[]
    l=size(thermal_ensemble)[1]

    ll=length(basis_int)
    T0=I(ll)
    Tpi=I(ll)

    @time begin
    T_power = I(ll) 
    # 计算 T 的幂并同时更新 T0 和 Tpi
    for i in 1:(N-1)
        T_power *= Trans  # 计算 T^i
        T0 += T_power
        Tpi += (-1)^i * T_power
    end
    end

    for i in 1:size(thermal_ensemble)[2]
        st=thermal_ensemble[:,i]
        if isapprox(st'*Trans*st, 1, atol=1e-6) 
            push!(total_thek0, st)
        elseif isapprox(st'*Trans*st, -1, atol=1e-6)
            push!(total_thekpi, st)
        else
            # stp,stm=vec2k0pi(N,st)
            stp=T0*st
            stp/=norm(stp)
            stm=Tpi*st
            stm/=norm(stm)
            push!(total_thek0, stp)
            push!(total_thekpi, stm)
        end
        
    end

    total_thek0 = hcat(total_thek0...)
    total_thekpi = hcat(total_thekpi...)

    return exact_scar, exact_scar_prime, total_thek0,total_thekpi
    
end
sep_scar_exact_translation(n::Int64, energy::Vector{Float64}, states::Matrix{Float64}) = sep_scar_exact_translation(BitStr{n, Int}, energy, states)