using Test
using PXPConstrained, BitBasis 
using ITensors
using Yao: arrayreg, density_matrix, von_neumann_entropy, expect, matblock
ITensors.set_warn_order(60)
include("../exm/FitEntEntScal.jl")

function translation_operator(n::Int, N::Int)
    dim = 2^N  # 希尔伯特空间的维度
    T = zeros(Float64, dim, dim)
    
    # 对每个基态应用平移
    for i in 0:(dim-1)
        # 将整数转换为二进制表示(位串)
        bits = digits(i, base=2, pad=N)
        
        # 执行循环平移
        shifted_bits = circshift(bits, n)
        
        # 将平移后的位串转换回整数
        j = sum(shifted_bits[k] * 2^(k-1) for k in 1:N)
        
        # 在适当位置放置矩阵元素(Julia是1-indexed)
        T[i+1, j+1] = 1.0
    end
    
    return T
end

function test_gene_scar(::Type{T}) where {N, T <: BitStr{N}}
    basis= PXP_basis(T)
    basis_int = [i.buf for i in basis]
    Trans = translation_operator(1,N)
    scar=gene_scar(N)
    scar1=storage(scar)
    scar1/=norm(scar1)
    scar2 = Trans*scar1

    # scar2 is the translation of scar1, scar2 is normalized
    @test norm(scar2) ≈ 1
   
    reg1 = arrayreg(scar1)
    reg2 = arrayreg(scar2)

    scark0=reg1+reg2
    scark0/=norm(scark0)
    scarkpi=reg1-reg2
    scarkpi/=norm(scarkpi)

    rhok0 = density_matrix(scark0, 1:div(N,2))
    rhokpi = density_matrix(scarkpi, 1:div(N,2))
    sk0=von_neumann_entropy(rhok0)
    skpi=von_neumann_entropy(rhokpi)
    c=expect(matblock(Trans), scark0)


    myscar1=scar1[basis_int.+1]
    myscar1/=norm(myscar1)
    myscar2=translation_matrix(T)*myscar1
    myscark0=myscar1+myscar2
    myscark0/=norm(myscark0)
    myscarkpi=myscar1-myscar2
    myscarkpi/=norm(myscarkpi)

    myrhok0=rdm_PXP(T, collect(1:div(N,2)), myscark0)
    myrhokpi=rdm_PXP(T, collect(1:div(N,2)), myscarkpi)
    mysk0=ee(myrhok0)
    myskpi=ee(myrhokpi)

    return c, sk0, mysk0, skpi, myskpi, myscar1, myscar2
end

@testset "scar separate" begin
    N=12
    T = BitStr{N, Int}
    proj=proj_FSA(T)
     # projector of subspace should be identity matrix
    @test isapprox(proj, I(size(proj)[1]), atol=1e-6)
   

    proj=proj_FSA2total(T)
    @test size(proj) == (322, 322)

    H = PXP_Ham(T)
    energy, states= eigen(H)
    scar, total_states = sep_scar_FSA(T, energy, states)

    splitlis = collect(1:N-1)
    EE_lis=zeros(length(splitlis))
    for m in eachindex(EE_lis)
        subrho=rdm_PXP(N, collect(1:splitlis[m]), scar)
        EE_lis[m]=ee(subrho)
    end

    cent, _= fitCCEntEntScal(EE_lis; mincut=1, pbc=true)
    @test isapprox(cent, 1.860152808765914, atol=0.01)

    Ob=randn((322, 322))
    proj=proj_Ob(energy, states, Ob)
    @test size(proj)==(22,22)

    scar, thermal = sep_scar_FSA(T, energy, states)
    total_st=hcat(scar, thermal)
    # The total_st should be orthonormal
    @test isapprox(total_st'*total_st, I(22))

    scar=gene_scar(N)
    scar1=storage(scar)

    basis= PXP_basis(T)
    basis_int = [i.buf for i in basis]
    scar1=scar1[basis_int.+1]
    @test isapprox(norm(H*scar1), 0, atol=1e-12) # This is the definition of a exact zero energy scar state.

    exact_scar,exact_scar_prime,thermal_ensemble=sep_scar_exact(T, energy, states)
    total_st = hcat(exact_scar, exact_scar_prime, thermal_ensemble)
    # The total_st should be orthonormal and idempotent
    @test isapprox(total_st'*total_st, I(24))
    P=total_st*total_st'
    @test isapprox(P^2, P)

    # cross-check the entropy of scar state
    c, sk0, mysk0, skpi, myskpi, myscar1, myscar2 = test_gene_scar(T)
    @test isapprox(sk0, mysk0)
    @test isapprox(skpi, myskpi)
    @test isapprox(c, 1)

    scar=gene_scar(6)
    scar1=storage(scar)

    basis= PXP_basis(BitStr{6,Int})
    basis_int = [i.buf for i in basis]
    scar1=scar1[basis_int.+1]
    scar2=translation_matrix(BitStr{6,Int})*scar1

    # for N<=6, scar1 and scar2 identical up to a overall coefficient.
    @test scar1 ≈ -scar2
end

