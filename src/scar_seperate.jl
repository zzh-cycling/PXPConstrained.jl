
function proj_FSA(N::Int64)
    """
    Generate the projection matrix for the FSA subspace, which should equals to identity matrix in FSA subspace. Meanwhile not equal to identity matrix in the total Hilbert space. Input N is the size of the system, return the projection matrix.
    """
    energy,states=eigen(PXP_FSA_Ham(N))

    return states*states'
end

function proj_FSAinTotal(N::Int64)
    # iso=iso_total2FSA(N)

    iso=load("/Users/cycling/Documents/projects/big_data/scar_thermal_FSA/iso_FSA/iso_total2FSA$(N).jld", "iso")

    myprint(stdout,"iso complete")
    energy,states=eigen(PXP_FSA_Ham(N))
    Proj=iso*states*states'*iso'
    return Proj
end

function sep_scar_FSA(N::Int64,energy::Vector{Float64},states::Matrix{Float64})
    indices=[index for (index,value) in enumerate(energy) if abs(value)<=1e-8]
    P_FSA=proj_FSAinTotal(N)

    myprint(stdout,"P_FSA complete")
    iso_0modes=states[:,indices]
    # Here need to note that for different basic MKL/openblas version, the degenrated subspace will generate different states. So we need to use the same kind of Package.
    # P_0modes=iso_0modes*iso_0modes'

    PPP_symmetrized = (iso_0modes'*P_FSA*iso_0modes + (iso_0modes'*P_FSA*iso_0modes)') / 2
    myprint(stdout,"PPP complete")
    vals, vecs = eigen(PPP_symmetrized)
    vecs=iso_0modes*vecs

    Inv=Inv_proj_matrix(N)
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
    # myprint(stdout,"PPP complete")
    # scar,thermal=vecs[:,end],vecs[:,1:end-1]

    # return scar,thermal
end

function sep_thermal(N,states,indices)
    P_FSA=proj_FSAinTotal(N)
    l=size(P_FSA,1)
    iso_0modes=states[:,indices]
    P_0modes=iso_0modes*iso_0modes'
    P_thermal=I(l)-P_FSA
    PPP_symmetrized = (P_0modes*P_thermal*P_0modes + (P_0modes*P_thermal*P_0modes)') / 2

    vals, thermalvecs = eigen(PPP_symmetrized)

    return  thermalvecs
end

function proj_Ob(energy::Vector{Float64},states::Matrix{Float64},Ob::Matrix{Float64})
    # Create the operator in the new basis of zero energy subspaces
    # Actually we find that the non correct overlap between the eigenstates and Z2states is due to the superposition of the degenrate scar states and thermal states in zero modes subspace. Not mainly due to the inversion symmetry.
    iso_0modes=states[:,indices]
    projected_matrix = iso_0modes'*Ob*iso_0modes
    return projected_matrix
end

function proj_Z2(energy::Vector{Float64}, states::Matrix{Float64})
    # Create the |Z2state><Z2state| in the new basis
    # Actually we find that the non correct overlap between the eigenstates and Z2states is due to the superposition of the degenrate scar states and thermal states in zero modes subspace. Not mainly due to the inversion symmetry.
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

function proj_Z2Zt2(N::Int64, states::Matrix{Float64},indices::Vector{Int64})
    # Create the operator in the new basis
    # Actually we find that the non correct overlap between the eigenstates and Z2states is due to the superposition of the degenrate scar states and thermal states in zero modes subspace. Not mainly due to the inversion symmetry.
    basis_int, basis_string = PXP_basis(N)
    z2tilde=findfirst(x->x==parse(Int, "01"^div(N,2), base=2),basis_int)+1
    projected_matrix = zeros(length(indices), length(indices))
    for i in eachindex(indices)
        for j in i: length(indices)
            u = real(states[:, indices[i]])
            v = real(states[:, indices[j]])
            # Because we are in the constrained subspace, the overlap between the states and Z2state is just the last element of the state. p[i,j]=u[end]*v[end] basically saying 
            #<state i|Z2state><Z2state|state j>
            
            projected_matrix[i, j] = u[end] * v[end]+u[z2tilde] * v[z2tilde]#if we choose u[1597] * v[1597] we get the same result, not only overlap and entanglement entropy.
            projected_matrix[j, i] = projected_matrix[i, j]
        end
    end
    return projected_matrix
end

function proj_invZ2(N::Int64, states::Matrix{Float64},indices::Vector{Int64})
    # Create the operator in the new basis
    # Actually we find that the non correct overlap between the eigenstates and Z2states is due to the superposition of the degenrate scar states and thermal states in zero modes subspace. Not mainly due to the inversion symmetry.
    basis_int, basis_string = PXP_basis(N)
    projected_matrix = zeros(length(indices), length(indices))
    for i in eachindex(indices)
        for j in i: length(indices)
            u = real(states[:, indices[i]])
            v = real(states[:, indices[j]])
            # Because we are in the constrained subspace, the overlap between the states and Z2state is just the last element of the state. p[i,j]=u[end]*v[end] basically saying 
            #<state i|Z2state><Z2state|state j>
            z2tilde=findfirst(x->x=="01"^div(N,2),basis_string)+1
            projected_matrix[i, j] = u[end] * v[end]+u[z2tilde] * v[z2tilde]+u[end] *v[z2tilde]+ u[z2tilde] *v[end]#if we choose u[1597] * v[1597] we get the same result, not only overlap and entanglement entropy.
            projected_matrix[j, i] = projected_matrix[i, j]
        end
    end
    return projected_matrix
end

function sep_scar_state(N::Int64,energy::Vector{Float64},states::Matrix{Float64})
    indices = [index for (index, value) in enumerate(energy) if abs(value) <= 1e-8]

    project_matrix = proj_Z2Zt2(N,states,indices);
    proj_val, proj_vec = eigen(project_matrix);
    scar1 = states[:,indices] * proj_vec[:,end];
    scar2= states[:,indices] * proj_vec[:,end-1];
    thermal_ensemble=states[:,indices] * proj_vec[:,1:end-2];

    return scar1,scar2,thermal_ensemble
end

function sep_scar_Ob(energy::Vector{Float64},states::Matrix{Float64},Ob::Matrix{Float64})
    indices = [index for (index, value) in enumerate(energy) if abs(value) <= 1e-8]

    project_matrix = proj_Ob(states,indices,Ob);
    proj_val, proj_vec = eigen(project_matrix);
    scar1 = states[:,indices] * proj_vec[:,end];
    scar2= states[:,indices] * proj_vec[:,end-1];
    thermal_ensemble=states[:,indices] * proj_vec[:,1:end-2];

    return scar1,scar2,thermal_ensemble
end


function superposition_scar(N::Int64, energy::Vector{Float64}, states::Matrix{Float64},lambda::Float64)
    scar, thermal_ensemble= sep_scar_FSA(N, energy,states)
    len=size(thermal_ensemble)[2]
    super_states=zeros(size(scar)[1],len)
    for i in 1:len 
        thermal=thermal_ensemble[:,i]
        state= (1-lambda) .*scar .+ lambda .* thermal
        if norm(state) > 0
            state= state/norm(state)
        end
        super_states[:,i]=state
    end
    
    return super_states
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

function vec2k0pi(N::Int64, state::Vector{Float64})
    basis_int, basis_string=PXP_basis(N)
    l=length(basis_int)
    T=translation_matrix(N)
    op1=I(l)
    op2=I(l)


    T_power = I(l)  
    # 计算 T 的幂并同时更新 op1 和 op2
    for i in 1:(N-1)
        T_power *= T  # 计算 T^i
        op1 += T_power
        op2 += (-1)^i * T_power
    end

    statek0=op1*state
    statek0/=norm(statek0)
    statekpi=op2*state
    statekpi/=norm(statekpi)
    return statek0,statekpi
end

function gram_schmidt(vectors::Matrix{Float64})
    n = size(vectors, 2)
    m = size(vectors, 1)

    orthogonal_vectors = zeros(Float64, m, n)

    for i in 1:n
        
        v_i = vectors[:, i]

        for j in 1:(i-1)
            proj_j = (dot(v_i, orthogonal_vectors[:, j]) / dot(orthogonal_vectors[:, j], orthogonal_vectors[:, j])) * orthogonal_vectors[:, j]
            v_i -= proj_j
        end

        v_i /= norm(v_i)
        orthogonal_vectors[:, i] = v_i
    end

    return orthogonal_vectors
end

function sep_scar_exact(N::Int64, energy::Vector{Float64}, states::Matrix{Float64})
    indices = [index for (index, value) in enumerate(energy) if abs(value) <= 1e-8]
    T = translation_matrix(N)
    scar=gene_scar(N)
    scar1=storage(scar)

    basis_int, basis_string= PXP_basis(N)
    scar1=scar1[basis_int.+1]
    scar1/=norm(scar1)

    scar2=T*scar1
    scar2/=norm(scar2)

    exact_scar= scar1+scar2
    exact_scar/=norm(exact_scar)

    exact_scar_prime=scar1-scar2
    exact_scar_prime/=norm(exact_scar_prime)

    thermal_ensemble=gram_schmidt(hcat(exact_scar,exact_scar_prime,states[:,indices]))[:,3:end]


    return exact_scar,exact_scar_prime,thermal_ensemble
    
end
