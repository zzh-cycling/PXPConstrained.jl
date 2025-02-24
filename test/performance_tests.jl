using BenchmarkTools
using Profile
using BitBasis, PXPConstrained

@profile rdm_PXP_K(BitStr{24, Int}, collect(1:12), ones(4341),0)
@time rdm_PXP_K(BitStr{26, Int}, collect(1:13), ones(10462),0)
Profile.print(format=:flat, mincount=40)
@time rdm_PXP_K(BitStr{28, Int}, collect(1:14), ones(25415),0)


@btime iso_total2FSA(BitStr{N, Int})
@btime iso_total2FSA(BitStr{N})
# Have difference

@btime PXP_Ham(BitStr{12, Int})

@btime PXP_Ham2(12)

@code_warntype PXP_Ham(BitStr{12, Int})
@code_warntype PXP_Ham2(12)
@btime PXP_basis(BitStr{16, Int})   
@btime PXP_basis2(16)
