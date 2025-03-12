This note is aimed to illustrate the issue of my rdm_K rdm_MSS function performance.

# Basics of dimension of matrix
## Constrained space
The dim of PXP model's Hilbert space after ruling out Rydberg blockade kinetically constrained states is a Fibonacci sequence

$$
\begin{cases}
&d_L^{OBC}=F_{L+2} \\
& d_L^{PBC}=F_{L-1}+F_{L+1}
\end{cases}
$$

Therefore, the dimension of the Hilbert space with L sites, $d_L$, satisfies the linear recurrence equation
$$
d_L = d_{L−1} + d_{L−2}
$$

The Hilbert space for L sites with PBC can be formed by taking the Hilbert space of the same number of sites with open boundary conditions and removing all configurations which both begin and end with •, hence, $d^{PBC}_L = d_L − d_{L−4}= F_{L−1} + F_{L+1}$.

## Symmetric basis

Then using the symmetry of this system to build symmetric basis out of constrained computational basis, for example take $k=0$ and inversion symmetric $I=1$ state as basis to represent the system.

1. using the representative states of each basis of one irreducible representation space to build the set of translational basis with different momentum.

$$
|\bar{n},k\rangle=\sum_{n=0}^{L-1} e^{\frac{i 2 \pi}{L}kn}T^{n}|\bar{n}\rangle
$$

For given momentum sector, we can write the basis out of representative basis, noticing \# representative states = # basis of $k=0$ sector.


2. Then we use inversion symmetry to simplify it further by reflection operator $R=R^{\dagger}$ defined as:

$$
RO_l R=O_{L-i-l}
$$

In the sector of $k=0, k=\pi$, $[T,R]=0$, because $RTR^{-1}=T^{-1}$ and $T=T^{-1}$, then we can continue divide it to $R=1$ sector and $R=-1$ sector.  

Finally we get $|S_{N=L/2}=O(L^{-1/2}2^L)|$, $|S_{N=L/2,k=0}=O(L^{-3/2}2^L)|$, are called the maximum symmetry sector (MSS)


## Matrix Diagonalization
Now we are intended to write the Hamiltonian $H_{k=0}$ and $H_{k=0,I=1}$. 

Then we find that the $k=0$ sector in Fibonacci sequence and inversion symmetry part is:

|size|Hilbert space dim|constrained Hilbert space dim|representatives/k=0|pure inversion|k=0,R=1 MSS|
|----|----|----|----|----|----|
|L=16 |65536|2207  |143 |55|99|
|L=24 |16777216|103682|4341  |377 |2359|
|L=26 |67108864|271443|10462 |610 |5536|
|L=28 |268435456|710647|25415 |987 |13201|
|L=30 |1073741824|1860498|62075 |1597 |  31836|
|L=32 |4294967296|4870847|152288 | 2584 |77436|

(# representatives+# pure inversion)/2= # MSS

and the number of translational invariant basis and the number of MSS differ when particle number is larger than 8

maximal fully diagonal system size is $L=32$, the biggest block of $L=34$ is 212534. As for the sparse matrix to calculate dynamics, we can scale up $L=38$.

For XXZ model, we could calculate the $L=24$.

# Calculate the Reduced density matrix

we want to calculate the reduced density matrix of a given pure state(in constrained space, further in symmetric constrained space.)
## Core formula
$$
\rho=\sum_k I \bigotimes \braket{k|\rho|k}\bigotimes I
$$

where $\ket{k}$ is the basis of traced out subsystem. The simplest code for total Hilbert space is:

```julia
# artificially divide system to left and right, trace out right.
n_right = n_total - n_left

# 计算子系统的维度
dim_left = 2^n_left
dim_right = 2^n_right
# 初始化约化密度矩阵
rho_reduced = zeros(ComplexF64, dim_left, dim_left)
# 遍历所有可能的子系统状态
for i in 1:dim_left
    for j in 1:dim_left
        # trace out common index k
        for k in 1:dim_right
            # 计算全局态的索引
            index_i = (i - 1) * dim_right + k
            index_j = (j - 1) * dim_right + k
            # 累加到约化密度矩阵
            rho_reduced[i, j] += psi[index_i] * conj(psi[index_j])
        end
    end
end
```

If we deal with constrained space, we can still follow the above procedure, just do the map from index of untraced out system to index of subsystem. We can have a example code(I do not quite understand all these lines):

```julia
reduced_dm = zeros(ET, (size, size))
cache = [Dict(zip(environ[idcs], idcs)) for idcs in sub_basis_index]
for i in 1:size, j in i:size
    val = zero(ET)
    for idr in sub_basis_index[i]
        r = environ[idr]
        target = get(cache[j], r, -1)
        if target != -1
            @inbounds val += state[idr]' * state[target]
        end
    end
    reduced_dm[i, j] = val
    reduced_dm[j, i] = reduced_dm[i, j]
```
which is the core part of the function of `rdm_PXP()` in `PXPConstrained`.

## Performance problem

The main problem we encounter is how to obtain the rdm of state from symmetric constrained space, for example, the states built on `PXP_K_basis()`,`PXP_MSS_basis()`. My method is find the mapping between K space and MSS space to original constrained space, then reuse the  `rdm_PXP()`. So the `rdm_PXP_K()` function contain two parts, one is `iso_total2K`, another is `rdm_PXP()`.

Here comes the issue [code](https://github.com/zzh-cycling/PXPConstrained/issues/3).

Theoretically, the method to solve this can still by idea given above, find the corresponding relation between index of K space and index of subsystem's K space. However, the problem lies that:

The total system has translational symmetry, meanwhile the subsystem artificially breaks such symmetry. We can not represent it in K space, but still in constrained space, the maximal total system size we deal with in merely $L=32$, the subsystem size can only be $L=16$, absolutely under our capacity. So, how can we directly find：

```mermaid
flowchart LR
A(K space) <--map relation--> B(Constrained subsystem space) 
```

instead of 

```mermaid
flowchart LR
A(K space) <--iso_total2K--> B(Constrained space) <--rdm_PXP--> C(Constrained subsystem space)
```

Where K space's state belike($L=4$):

$\ket{1,k=0}=\frac{1}{2}(\ket{1}+\ket{2}+\ket{4}+\ket{8})$

From the perspective of algorithm, we may only consider the representative state $\ket{1}$'s map to subsystem, them iterate them by translation?
