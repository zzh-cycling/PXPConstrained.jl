using BenchmarkTools
using Profile
using BitBasis, PXPConstrained

@profile rdm_PXP_K(BitStr{24, Int}, collect(1:12), ones(4341),0)
@time rdm_PXP_K(BitStr{26, Int}, collect(1:13), ones(10462),0)
Profile.print(format=:flat, mincount=40)
@time rdm_PXP_K(BitStr{28, Int}, collect(1:14), ones(25415),0)
