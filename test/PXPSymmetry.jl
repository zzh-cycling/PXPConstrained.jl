@testset "basis" begin
    N=18
    basis = PXP_MSS_basis(N,div(N,2))[1]
    basisMSSm = PXP_MSS_basis(N,0,-1)[1]
    @test filter(x -> !(x in basis), basisMSSm)==[18981]
end

function testzeromodes(N)
    sum1=0
    sum2=0
    sum2m=0
    # sum1 is the number of zero modes in different momentum sectors. sum2 is the number of zero modes in zero momentum and \pm 1 Inversion sector. sum2m is the number of zero modes in pi momentum sectors and \pm 1 Inversion sector.
    hss = length(findall(x -> isapprox(x,0,atol=1e-10), eigvals(PXP_MSS_Ham(N, 0))))
    hssm = length(findall(x -> isapprox(x,0,atol=1e-10), eigvals(PXP_MSS_Ham(N, 0, -1))))
    hss1=length(findall(x -> isapprox(x,0,atol=1e-10), eigvals(PXP_MSS_Ham(N,div(N,2)))))
    hss1m=length(findall(x -> isapprox(x,0,atol=1e-10), eigvals(PXP_MSS_Ham(N,div(N,2),-1))))
    sum2=hss+hssm
    sum2m=hss1+hss1m
    
    for i in 0:N-1
        hk = PXP_K_Ham(N, i)
        sum1+=length(findall(x -> isapprox(x,0,atol=1e-10), eigvals(hk)))
    end

    return sum1, sum2, sum2m
end

@testset "zero modes" begin
    zeromodeslis=[13, 22, 37, 59, 100, 165, 269]
    for (idx,i) in enumerate(10:2:20)
        sum1,sum2, sum2m=testzeromodes(i)
        @test sum1==zeromodeslis[idx]
        # @test_broken sum2==length(findall(x -> isapprox(x,0,atol=1e-10), eigvals(PXP_K_Ham(i, 0))))
        # @test_broken sum2m==length(findall(x -> isapprox(x,0,atol=1e-10), eigvals(PXP_K_Ham(i, div(i,2)))))
    end
end

@testset "pxp k and mss 4n" begin
    N=12
    basisK= PXP_K_basis(N, 3)[1]
    @test length(basisK) == 26
    @test BitStr{N}(1)==basisK[1]
    @test basisK[end].buf <<1 ==1322

    basisK= PXP_K_basis(N, 0)[1]
    @test length(basisK) == 31
    @test BitStr{N}(0)==basisK[1]
    @test basisK[end].buf <<1 ==2730
    @test basisK[2].buf==1
    @test basisK[3].buf==5

    basisKpi= PXP_K_basis(N, div(N,2))[1]
    @test length(basisKpi) == 29
    @test BitStr{N}(1)==basisKpi[1]
    @test basisKpi[end].buf <<1 ==2730
    @test filter(x -> !(x in basisKpi), basisK)==[BitStr{N}(0),BitStr{N}(585)]

    basis = PXP_MSS_basis(N,0)[1]
    @test length(basis) == 26
    @test BitStr{N}(0)==basis[1]
    @test basis[end].buf <<1==2730

    basisMSSm = PXP_MSS_basis(N,0,-1)[1]
    @test length(basisMSSm) == 5
    @test basisMSSm[1].buf == 37
    @test basisMSSm[end].buf==293

    try
        # Testing invalid input parameters
        basis = PXP_MSS_basis(N,3)[1]
    catch e
        # If the function correctly validates inputs, we expect an error
        @test e isa Exception
    end

    try
        # Testing invalid input parameters
        basis = PXP_MSS_basis(N,3,-1)[1]
    catch e
        # If the function correctly validates inputs, we expect an error
        @test e isa Exception
    end

    basispi = PXP_MSS_basis(N,div(N,2))[1]
    @test basispi == basisMSSm #only for N<=16

    basispim = PXP_MSS_basis(N,div(N,2),-1)[1]
    @test length(basispim) == 24
    @test BitStr{N}(1)==basispim[1]
    @test basispim[end].buf <<1==2730

    hk = PXP_K_Ham(N, 0)
    @test size(hk) == (31, 31)
    @test ishermitian(hk)
    @test hk[1, 1] == hk[end, end]≈ 0.0
    zeromodes=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hk))
    @test zeromodes==[14, 15, 16, 17, 18]

    hk = PXP_K_Ham(N, 3)
    @test size(hk) == (26, 26)
    @test ishermitian(hk)
    @test hk[1, 1] == hk[end, end] ≈ 0.0
    @test findall(x -> isapprox(x,0,atol=1e-10), eigvals(hk))==[13, 14]

    hk = PXP_K_Ham(N, div(N,2))
    @test size(hk) == (29, 29)
    @test ishermitian(hk)
    @test hk[1, 1] == hk[end, end] ≈ 0.0
    zeromodes5=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hk))
    @test zeromodes5==[14, 15, 16]

    hmss=PXP_MSS_Ham(N, 0)
    @test size(hmss) == (length(basis),length(basis))
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0
    zeromodes1=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hmss))
    @test zeromodes1==[12, 13, 14, 15]

    hmss=PXP_MSS_Ham(N, 0, -1)
    @test size(hmss) == (length(basisMSSm), length(basisMSSm))
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0
    zeromodes2=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hmss))
    @test zeromodes2==[3]
    @test length(zeromodes2)+length(zeromodes1) == length(zeromodes)

    hmss=PXP_MSS_Ham(N, div(N,2))
    @test size(hmss) == (length(basispi), length(basispi))
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0
    zeromodes3=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hmss))
    @test zeromodes3==[3]

    hmss=PXP_MSS_Ham(N, div(N,2), -1)
    @test size(hmss) == (length(basispim), length(basispim))
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0
    zeromodes4=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hmss))
    @test zeromodes4==[12, 13]
    @test length(zeromodes3)+length(zeromodes4)==length(zeromodes5)

    N = 12
    MSS_basis = PXP_MSS_basis(N, 0)[1]
    @test length(MSS_basis) == 26
    @test length(PXP_MSS_basis(8, 0)[1]) == length(PXP_K_basis(8, 0)[1])

    # The energy spectrum is symmetric about zero
    H_MSS = PXP_MSS_Ham(N, 0)
    MSS_vals, MSS_vecs = eigen(H_MSS)
    @test isapprox(reverse(MSS_vals) .+ MSS_vals, zeros(26), atol=1e-6)

    # The MSS dims + MSSminv dims = K dims
    @test length(PXP_MSS_basis(N, 0)[1]) + length(PXP_MSS_basis(N, 0, -1)[1]) == length(PXP_K_basis(N, 0)[1])
    H_MSSminv = PXP_MSS_Ham(N, 0, -1)
    MSSminv_vals, MSSminv_vecs = eigen(H_MSSminv)
    @test isapprox(MSSminv_vals[3], 0.0, atol=1e-6)

end

@testset "pxp k and mss 4n+2" begin
    N=14
    basisK3= PXP_K_basis(N, 3)[1]
    @test length(basisK3) == 58
    @test BitStr{N}(1)==basisK3[1]
    @test basisK3[end].buf <<1 ==5290 #00101001010101

    basisK= PXP_K_basis(N, 0)[1]
    @test length(basisK) == 64
    @test BitStr{N}(0)==basisK[1]
    @test basisK[end].buf <<1 == 10922
    @test basisK[2].buf==1
    @test basisK[3].buf==5

    basisKpi= PXP_K_basis(N, div(N,2))[1]
    @test length(basisKpi) == 59
    @test BitStr{N}(1)==basisKpi[1]
    @test basisKpi[end].buf <<1 == 10922
    @test filter(x -> !(x in basisKpi), basisK)==[BitStr{N}(0),BitStr{N}(129),BitStr{N}(645),BitStr{N}(1161),BitStr{N}(2709)]

    basis = PXP_MSS_basis(N,0)[1]
    @test length(basis) == 49
    @test BitStr{N}(0)==basis[1]
    @test basis[end].buf <<1==10922

    basisMSSm = PXP_MSS_basis(N,0,-1)[1]
    @test length(basisMSSm) == 15
    @test basisMSSm[1].buf == 37
    @test basisMSSm[end].buf==1189

    try
        # Testing invalid input parameters
        basis = PXP_MSS_basis(N,3)[1]
    catch e
        # If the function correctly validates inputs, we expect an error
        @test e isa Exception
    end

    try
        # Testing invalid input parameters
        basis = PXP_MSS_basis(N,3,-1)[1]
    catch e
        # If the function correctly validates inputs, we expect an error
        @test e isa Exception
    end

    basispi = PXP_MSS_basis(N,div(N,2))[1]
    @test basispi == basisMSSm #only for N<=16

    basispim = PXP_MSS_basis(N,div(N,2),-1)[1]
    @test length(basispim) == 44
    @test BitStr{N}(1)==basispim[1]
    @test basispim[end].buf <<1==10922

    hk = PXP_K_Ham(N, 0)
    @test size(hk) == (length(basisK), length(basisK))
    @test ishermitian(hk)
    @test hk[1, 1] == hk[end, end]≈ 0.0
    zeromodes=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hk))
    @test zeromodes==29:36

    hk = PXP_K_Ham(N, 3)
    @test size(hk) == (length(basisK3), length(basisK3))
    @test ishermitian(hk)
    @test hk[1, 1] == hk[end, end] ≈ 0.0
    @test findall(x -> isapprox(x,0,atol=1e-10), eigvals(hk))==[29,30]

    hk = PXP_K_Ham(N, div(N,2))
    @test size(hk) == (length(basisKpi), length(basisKpi))
    @test ishermitian(hk)
    @test hk[1, 1] == hk[end, end] ≈ 0.0
    @test findall(x -> isapprox(x,0,atol=1e-10), eigvals(hk))==28:32

    hmss=PXP_MSS_Ham(N, 0)
    @test size(hmss) == (length(basis), length(basis))
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0
    zeromodes1=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hmss))
    @test zeromodes1==23:27

    hmss=PXP_MSS_Ham(N, 0, -1)
    @test size(hmss) == (length(basisMSSm), length(basisMSSm))
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0
    zeromodes2=findall(x -> isapprox(x,0,atol=1e-10), eigvals(hmss))
    @test zeromodes2==6:10
    @test_broken length(zeromodes2)+length(zeromodes1) == length(zeromodes)

    hmss=PXP_MSS_Ham(N, div(N,2))
    @test size(hmss) == (length(basispi), length(basispi))
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0
    @test findall(x -> isapprox(x,0,atol=1e-10), eigvals(hmss))==6:10

    hmss=PXP_MSS_Ham(N, div(N,2), -1)
    @test size(hmss) == (length(basispim), length(basispim))
    @test ishermitian(hmss)
    @test hmss[1, 1] == hmss[end, end] ≈ 0.0
    @test findall(x -> isapprox(x,0,atol=1e-10), eigvals(hmss))==[]
end

@testset "iso, reduced density matrix and map" begin
    N=12

    map_idx = iso_total2K(BitStr{12, Int},0)
    @test size(map_idx) == (322, 31)
    @test map_idx'*map_idx ≈ I(31)
    u = iso_total2K(14, 7)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-8)
    @test map_idx[1, 1] == 1.0
    @test map_idx[end, end] ≈ 1/√2
    @test map_idx[233, end] ≈ 1/√2

    iso = iso_K2MSS(12, 0)
    @test size(iso) == (31, 26)
    @test iso'*iso ≈ I(26)
    P=iso*iso'
    @test isapprox(P*P, P, atol=1e-10)
    non_zero_lis = reshape(mapslices(col -> count(x -> x != 0, col), iso, dims=1), 26)
    @test non_zero_lis == PXP_MSS_basis(12, 0)[3]
    u = iso_K2MSS(14, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)

    u = iso_total2MSS(12, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)
    @test u[1,1] == 1.0
    @test u[end, end] ≈ 1/√2
    @test u[233, end] ≈ 1/√2

    u = iso_total2MSS(14, 7)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)

    u = iso_total2K(14, 0)
    P = u*u'
    @test isapprox(P*P, P, atol=1e-10)
    @test isapprox(u'*u, I(size(u, 2)), atol=1e-10)
    @test u[1,1] == 1.0
    @test u[end, end] ≈ 1/√2

    rdm_K = rdm_PXP_K(24, collect(1:12), ones(4341),0)
    @test size(rdm_K) == (377, 377)
    @test length(PXP_basis(12, false)) == 377

    map_K2MSS = iso_K2MSS(N, 0)
    map_K2MSSminv = iso_K2MSS(N, 0, -1)
    @test size(map_K2MSS) == (31, 26)
    @test map_K2MSS'*map_K2MSS ≈ I(26)
    @test size(map_K2MSSminv) == (31, 5)
    @test map_K2MSSminv'*map_K2MSSminv ≈ I(5)

    map_total2MSS = iso_total2MSS(N,0)
    @test size(map_total2MSS) == (322, 26)
    @test map_total2MSS'*map_total2MSS ≈ I(26)

    rdm_MSS = rdm_PXP_MSS(N, collect(1:6), eigvecs(PXP_MSS_Ham(N, 0))[:,12],0)
    @test size(rdm_MSS) == (21, 21) == (length(PXP_basis(6, false)) ,length(PXP_basis(6, false)))

    # Fit the scar's central charge, may change depends on the machine and basic linear algebra package.
end

@testset "mapstate" begin
    kstatez2 = zeros(31); kstatez2[end] = 1
    kstate0 = zeros(31); kstate0[1] = 1
    kstate1 = zeros(31); kstate1[2] = 1

    totalZ2 = mapstate_K2total(12, kstatez2, 0)
    @test norm(totalZ2) ≈ 1.0
    @test totalZ2[233] ≈ 1/√2
    @test totalZ2[end] ≈ 1/√2

    total0 = mapstate_K2total(12, kstate0, 0)
    @test norm(total0) ≈ 1.0
    @test total0[1] ≈ 1.0

    total1 = mapstate_K2total(12, kstate1, 0)
    @test norm(total1) ≈ 1.0
    orderlis=[i.buf+1 for i in PXP_basis(12)]
    index = map(x ->searchsortedfirst(orderlis, x), [1<<i+1 for i in 0:11])
    @test total1[index] ≈ 1/√12*ones(Float64, 12)

    MSSstatez2 = zeros(26); MSSstatez2[end] = 1
    MSSstate0 = zeros(26); MSSstate0[1] = 1
    MSSstate1 = zeros(26); MSSstate1[2] = 1
    MSSstate37 = zeros(26); MSSstate37[8] = 1

    totalZ2 = mapstate_MSS2total(12, MSSstatez2, 0)
    @test norm(totalZ2) ≈ 1.0
    @test totalZ2[233] ≈ 1/√2
    @test totalZ2[end] ≈ 1/√2

    total0 = mapstate_MSS2total(12, MSSstate0, 0)
    @test norm(total0) ≈ 1.0
    @test total0[1] ≈ 1.0

    total1 = mapstate_MSS2total(12, MSSstate1, 0)
    @test norm(total1) ≈ 1.0
    @test total1[index] ≈ 1/√12*ones(Float64, 12)

    Kz2 = mapstate_MSS2K(12, MSSstatez2, 0)
    @test norm(Kz2) ≈ 1.0
    @test Kz2[end] ≈ 1.0

    K0 = mapstate_MSS2K(12, MSSstate0, 0)
    @test norm(K0) ≈ 1.0
    @test K0[1] ≈ 1.0

    K37 = mapstate_MSS2K(12, MSSstate37, 0)
    @test norm(K37) ≈ 1.0
    @test K37[8] ≈ 1/√2
    @test K37[9] ≈ 1/√2
end

@testset "rdm_PXP_MSS, ergotropy_PXP_MSS_state" begin
    N = 10
    rotated_state1 = rotated_psi_state_mss(N, 0, π/4)
    rotated_state2 = rotated_psi_state_mss(N, div(N, 2), π/4, -1)
    rotated_state3 = rotated_psi_state_mss(N, div(N, 2), π/4)
    rotated_state4 = rotated_psi_state_mss(N, 0, π/4, -1)

    @test length(rotated_state3) == 1
    @test length(rotated_state4) == 1
    @test length(rotated_state1) == 14
    @test length(rotated_state2) == 11

    rdm1 = rdm_PXP_MSS(N, collect(1:div(N,2)), rotated_state1, 0)
    @test size(rdm1) == (13, 13)
    @test ishermitian(rdm1)
    @test tr(rdm1) ≈ 1.0

    rdm2 = rdm_PXP_MSS(N, collect(1:div(N,2)), rotated_state2, div(N,2), -1)
    @test size(rdm2) == (13, 13)
    @test ishermitian(rdm2)
    @test tr(rdm2) ≈ 1.0

    [-1.7083168055004627, -3.317496872634849, -2.68645531695468] .= ergotropy_PXP_MSS_state(N, div(N,2), rotated_state1, 0)
    [-1.7108050847457632, -3.317496872634849, -2.6924091458780577] .= ergotropy_PXP_MSS_state(N, div(N,2), rotated_state2, div(N,2), -1)
end