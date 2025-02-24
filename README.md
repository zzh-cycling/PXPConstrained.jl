# PXPConstrained

In constrained Hilbert considering the symmetry to do Exact Diagonalization. 

## Introduction

PXPConstrained is a a powerful exact numerical calculation package using Exact Diagonalization to deal with constrained quantum many-body problems, and leveraging the high performance of Julia. Especially, we use PXP model as example. It provides following functions:

- Generate (symmetric) basis in constrained Hilbert space, including Inversion symmetry, Translation symmetry. 
- Provide different physical quantity in constrained Hilbert space, such as `EE`(Entanglement Entropy), `TMI`(Tri-paritie Mutual Information), `MI`(Mutual Information), `QFI`(Quantum Fisher information).
- Separate the hybridized denegerated scar and thermal states in zero energy subspace.

## Try your first program

```julia
    N=12
    state=BitStr{N, Int}(0)
    basis = PXP_basis(BitStr{N, Int})
    H=PXP_Ham(BitStr{N, Int})
    rdm = rdm_PXP(BitStr{N, Int}, collect(1:6), state)
```

## Installation

```julia
pkg> add PXPConstrained
```

## Documentation

None