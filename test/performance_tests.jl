using BenchmarkTools
using Profile
using BitBasis, PXPConstrained

@profile rdm_PXP_K(BitStr{24, Int}, collect(1:12), ones(4341),0) # 50ms 
Profile.print(format=:flat, mincount=40)

@btime rdm_PXP(24, collect(1:12), ones(103682)) # 29ms
@btime rdm_PXP(BitStr{26, Int}, collect(1:13), ones(271443)) # 117ms
@btime rdm_PXP(BitStr{28, Int}, collect(1:14), ones(710647)) # 466ms
@btime rdm_PXP(30, collect(1:15), ones(1860498)) # 1.9s

@btime iso_total2K(24,0) # 357ms
@btime rdm_PXP_K(BitStr{24, Int}, collect(1:12), ones(4341),0) # 626ms
@btime rdm_PXP_K(BitStr{26, Int}, collect(1:13), ones(10462),0) #5.4s
@btime rdm_PXP_K(BitStr{28, Int}, collect(1:14), ones(25415),0)

@btime mapstate_total2K(BitStr{24, Int}, 0, ones(4341)) # 1.1ms

@btime iso_total2FSA(BitStr{N, Int})
@btime iso_total2FSA(BitStr{N})
# Have difference

@btime PXP_Ham(BitStr{12, Int})

@btime PXP_Ham2(12)

@code_warntype PXP_Ham(BitStr{12, Int})
@code_warntype PXP_Ham2(12)
@btime PXP_basis(BitStr{16, Int})   
@btime PXP_basis2(16)
