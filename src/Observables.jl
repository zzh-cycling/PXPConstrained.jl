include("PXP_basis.jl")
include("ptrace.jl")


function EE(subrm::Matrix{Float64})
    #  subrm=qi.ptrace(state*state',[2 for i in 1:N],[i for i in l+1:N])
     spectrum=eigvals(subrm)
     EE=0
     for i in eachindex(spectrum)
         v=abs(spectrum[i])
             if v>1e-8
                 EE+=-v*log(v)
             end
     end

     return EE
end

function EE_scaling_fig(N::Int64, state::Vector{Float64})
    splitlis=Vector(1:N-1)
    EElis=EE_PXP_state(N,splitlis,state)
    cent, fig=fitCCEntEntScal(EElis; mincut=1, pbc=true)

    return cent, fig
end

function EE_PXP_idx(N::Int64,splitlis::Vector{Int64}, idx::Int64)
    """
    Calculate the EE list for a given index of scar state
    """
    energy, states= eigen(PXP_Ham(N))
    idx_state=states[:,idx]
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subidx=rdm_PXP(N,collect(1:splitlis[m]), idx_state)
        EE_lis[m]=EE(subidx)
    end
    return EE_lis
end

function EE_PXP_idx_simplified(N::Int64,splitlis::Vector{Int64}, idx)
    """
    Simplified version of EE_PXP_idx, only calculate half the EE list
    """
    energy, states= eigen(PXP_Ham(N))
    idx_state=states[:,idx]
    EE_lis=zeros(div(length(splitlis)+1,2))
    for m in 1:div(length(splitlis)+1,2)
        subidx=rdm_PXP(N,collect(1:splitlis[m]), idx_state)
        EE_lis[m]=EE(subidx)
    end
    EE_lis=[EE_lis; sort(EE_lis[1:div(length(splitlis)-1,2)],rev=true)]
    return EE_lis
end

function EE_PXP_state(N::Int64,splitlis::Vector{Int64},state::Vector{Float64})
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subscar=rdm_PXP(N,collect(1:splitlis[m]), state)
        EE_lis[m]=EE(subscar)
    end
    return EE_lis
end

function Inv_proj_matrix(N::Int64, inv::Int64=1)
    basis_int, basis_string=PXP_basis(N)
    l=length(basis_int)
    Imatrix=zeros((l,l))
    reversed_basis_string=inversion_strlis(basis_string)
    for i in eachindex(basis_int)
        output=reversed_basis_string[i]
        j=findfirst(x -> x==parse(Int, output, base=2),basis_int)
        Imatrix[i,j]+=1.0
    end
    # positive_inversion_states=Inv_proj_matrix(16,1)*states[:,1075]
    # negative_inversion_states=Inv_proj_matrix(16,-1)*states[:,1075]
    # if inv==1
    #     return (I(l)+Imatrix)/2
    # end
    # if inv==-1
    #     return (I(l)-Imatrix)/2
    # end
    return Imatrix
end

function Mutual_information(N::Int64, state::Vector{Float64}, subsystems::Tuple{Vector{Int64}, Vector{Int64}})
    A, B = subsystems
    # MI formula defined as: I(A:B) = S_A + S_B - S_AB
    # Calculate the reduced density matrices
    # ρ_A = ptrace(N, A, state)
    # ρ_B = ptrace(N, B, state)
    # ρ_AB = ptrace(N, vcat(A,B), state)
    ρ_A = rdm_PXP(N, A, state)
    ρ_B = rdm_PXP(N, B, state)
    ρ_AB = rdm_PXP(N, vcat(A,B), state)
    # Calculate the Von Neumann entropies
    S_A = EE(ρ_A)
    S_B = EE(ρ_B)
    S_AB = EE(ρ_AB)
    # Calculate the mutual information
    I_AB = S_A + S_B - S_AB
    return I_AB
    
end

function Tri_mutual_information(N::Int64, state::Vector{Float64}, subsystems::Tuple{Vector{Int64}, Vector{Int64}, Vector{Int64}})
    A, B, C = subsystems
    # TMI formula defined as: I(A:B:C) = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    # Calculate the reduced density matrices
    # ρ_A = ptrace(N, A, state)
    # println(ρ_A==ρ_A',tr(ρ_A))
    # ρ_B = ptrace(N, B, state)
    # println(ρ_B==ρ_B',tr(ρ_B))
    # ρ_C = ptrace(N, C, state)
    # println(ρ_C==ρ_C',tr(ρ_C))

    # ρ_AB = ptrace(N, vcat(A,B), state)
    # println(ρ_AB==ρ_AB',tr(ρ_AB))
    # ρ_BC = ptrace(N, vcat(B,C), state)
    # println(ρ_BC==ρ_BC',tr(ρ_BC))
    # ρ_AC = ptrace(N, vcat(A,C), state)
    # println(ρ_AC==ρ_AC',tr(ρ_A))

    # ρ_ABC = rdm_PXP(N, vcat(A,B,C), state)
    # println(ρ_ABC==ρ_ABC',tr(ρ_ABC))
    # @time begin 
    ρ_A = rdm_PXP(N, A, state)
    # myprint(stdout,"rho_A complete")
    ρ_B = rdm_PXP(N, B, state)
    # myprint(stdout,"rho_B complete")
    ρ_C = rdm_PXP(N, C, state)
    # myprint(stdout,"rho_C complete")

    ρ_AB = rdm_PXP(N, vcat(A,B), state)
    # myprint(stdout,"rho_AB complete")
    ρ_BC = rdm_PXP(N, vcat(B,C), state)
    # myprint(stdout,"rho_BC complete")
    ρ_AC = rdm_PXP(N, vcat(A,C), state)
    # myprint(stdout,"rho_AC complete")

    ρ_ABC = rdm_PXP(N, vcat(A,B,C), state)
    # end
    # myprint(stdout,"rho complete")
    
    # Calculate the Von Neumann entropies
    # @time
    S_A = EE(ρ_A)
    # myprint(stdout,"SA EE complete")
    S_B = EE(ρ_B)
    S_C = EE(ρ_C)
    S_AB = EE(ρ_AB)
    # myprint(stdout,"SAB EE complete")
    S_BC = EE(ρ_BC)
    S_AC = EE(ρ_AC)
    S_ABC = EE(ρ_ABC)
    # end
    # myprint(stdout,"EE complete")

    # Calculate the mutual information
    I_ABC = S_A + S_B + S_C - S_AB - S_BC - S_AC + S_ABC
    
    return I_ABC
end

function QFI_dwd(N::Int64, Ob::Vector{Float64}, state::Vector{Float64})
    
    DeltaOb=state'*(Ob.^2 .*state)-(state'*(Ob.*state))^2
    # Calculate the Quantum Fisher Information
    # For spin 1/2, w/o 4
    F_Q = DeltaOb

    return F_Q
end

function Quantumfisherinfo(N::Int64, Ob::Matrix{Float64}, state::Vector{Float64})
    rho=state*state'
    DeltaOb=tr(rho*Ob^2)-tr(rho*Ob)^2
    # Calculate the Quantum Fisher Information
    # For spin 1/2, w/o 4
    F_Q = DeltaOb

    return F_Q
end

function domain_wall(N::Int64,pbc::Bool=true)
    """
    :param N: Number of sites
    :return: domain_wall_density operator
    The eigenvectors of this operator are going from -N to N, increasing by 2, totally N+1 eigenvectors. Number of each eigenvalues is N choose k, where k is the number of domain walls when we consider total Hilbert space. Defined as sum_i Z_i = (-1)^(i+1) * Z_i.
    """
    
    basis_int, basis_string = PXP_basis(N,pbc)
    l=length(basis_int)
    Dwd = zeros((l, l))

    function Z_i(str::String)
        return [c == '0' ? 1 : -1 for c in str]
    end

    function sum_Z_i(str::String)
        Z = Z_i(str)
        return sum((-1)^(i+1) * Z[i] for i in 1:length(Z))
    end

    for (idx, str) in enumerate(basis_string)
        Dwd[idx, idx] = sum_Z_i(str)
    end

    return Dwd
end

function domain_wall_density(N::Int64,pbc::Bool=true)
    """
    :param N: Number of sites
    :return: domain_wall_density operator
    The eigenvectors of this operator are going from -N to N, increasing by 2, totally N+1 eigenvectors. Number of each eigenvalues is N choose k, where k is the number of domain walls when we consider total Hilbert space.
    """
    
    basis_int, basis_string = PXP_basis(N,pbc)
    l=length(basis_int)
    Dwd = zeros((l, l))

    function Z_i(str::String)
        return [c == '0' ? 1 : -1 for c in str]
    end

    function sum_Z_i(str::String)
        Z = Z_i(str)
        return sum((-1)^(i+1) * Z[i] for i in 1:length(Z))
    end

    for (idx, str) in enumerate(basis_string)
        Dwd[idx, idx] = sum_Z_i(str)./N
    end

    return Dwd
end

function particlenumber(N::Int64,pbc::Bool=true)
    """
    :param N: Number of sites
    :return: Particle number operator
    """
    basis_int, basis_string = PXP_basis(N,pbc)
    l=length(basis_int)
    P = zeros((l, l))
    for (idx, str) in enumerate(basis_string)
        P[idx, idx] = count_ones(parse(Int, str, base=2))
    end

    return P
    
end

function on_siten(N::Int64, i::Int64,pbc::Bool=true)
    """
    :param N: Number of sites
    :return: Particle number operator
    """
    basis_int, basis_string = PXP_basis(N,pbc)
    l=length(basis_int)
    P = zeros((l, l))
    for (idx, str) in enumerate(basis_string)
        P[idx, idx] += (parse(Int, str[i], base=2))
    end

    return P
    
end
on_siten(4,2)
energy, states=eigen(PXP_Ham(16))

scar_indexlis16=[1, 2, 9, 27, 82, 202, 408, 728, 1075, 1480, 1800, 2006, 2126, 2181, 2199, 2206, 2207];
for j in scar_indexlis16
    sum=0
    st=states[:,j]
    for i in 1:16
        sum+=st'*on_siten(16,i)*st
    end
    println(energy[j], ",",sum/=16)
end


st=states[:,1100];
energy[1075]
sum=0
for i in 1:16
    sum+=st'*on_siten(16,i)*st
end
sum/=16
function swap_bits(bitstring::String)
    # 将比特字符串转换为字符数组
    bits = collect(bitstring)
    # 获取长度
    n = length(bits)
    
    # 创建一个新的数组来存储调换后的结果
    swapped_bits = copy(bits)
    
    # 遍历每两个元素进行调换
    for i in 1:2:n-1
        swapped_bits[i], swapped_bits[i+1] = swapped_bits[i+1], swapped_bits[i]
    end
    
    # 将字符数组转换回字符串并返回
    return String(swapped_bits)
end

function chiral(N::Int64, pbc::Bool=true)
    """
    :param N: Number of sites
    :return: Chiral operator
    N=6 found no index.
    """
    basis_int, basis_string = PXP_basis(N,pbc)
    l=length(basis_int)
    C = zeros((l, l))
    for (idx, str) in enumerate(basis_string)
        m=parse(Int, swap_bits(str), base=2)
        index=findfirst(x->x==m,basis_int)
        C[idx, index] +=1
    end

    return C
end