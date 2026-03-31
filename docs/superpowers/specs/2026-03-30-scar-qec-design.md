# Quantum Error Correction via Quantum Many-Body Scars

## Overview

Design a quantum error correcting code where the code space is the scar subspace of the PXP model, leveraging the algebraic structure (approximate SU(2) from FSA) for passive error protection. This work is **primarily based on** the ETH-AQECC framework (Brandão et al., arXiv:1710.04631), with methodology from CFT codes (Sang & Zou, arXiv:2406.09555).

## Core Insight

**QMBS as "anti-ETH" codes:** The Brandão et al. paper shows ETH eigenstates form AQECC because they are locally indistinguishable. Scar states violate ETH but may still form AQECC via a *different mechanism*: the approximate $\mathfrak{su}(2)$ algebraic structure from FSA, analogous to how CFT codes use conformal symmetry.

## Motivation

- **ETH codes** (Brandão et al.): Chaotic eigenstates form AQECC with $k = \Omega(N)$, $\varepsilon \sim e^{-cN}$
- **CFT codes** (Sang-Zou): Low-energy subspace with conformal symmetry; threshold $\Delta_{\min} > 1/2$
- **Scar codes** (this work): Use scar subspace with approximate SU(2) from FSA — *neither* ETH nor low-energy
- **Key question:** Do scars inherit AQECC properties from their algebraic structure despite violating ETH?

## Comparison of Three Frameworks

| Property | ETH Code (Brandão et al.) | CFT Code (Sang-Zou) | Scar Code (This Work) |
|----------|---------------------------|---------------------|----------------------|
| **Code space** | Finite energy density eigenstates | Low-energy eigenstates | Scar subspace $\mathcal{H}_{\rm scar}$ |
| **Physical system** | Chaotic many-body systems | Critical spin chains (e.g., TFIM) | PXP model (Rydberg arrays) |
| **Protection mechanism** | Local indistinguishability (ETH) | Conformal symmetry, scaling dimensions | Approximate $\mathfrak{su}(2)$ from FSA |
| **Thermalization** | Satisfies ETH | N/A (ground states) | **Violates ETH** (key distinction) |
| **Code rate $k$** | $\Omega(N)$ | $\Omega(\log\log n)$ | TBD (likely $O(1)$ to $O(\log L)$) |
| **Code distance $d$** | Determined by ETH locality | $\propto$ scaling dimension | TBD from FSA structure |
| **Error $\varepsilon$** | $\sim e^{-cN}$ | $\sim n^{-(2\Delta-1)}$ | TBD |
| **Threshold condition** | Always (for ETH systems) | $\Delta_{\min} > 1/2$ | TBD (FSA approximation quality?) |
| **Knill-Laflamme** | Satisfied via ETH diagonal/off-diagonal | Approximate via scaling | **To be verified** |
| **Energy location** | Middle of spectrum | Bottom of spectrum | Middle of spectrum (like ETH) |
| **Entanglement** | Volume law | Area law (ground state) | **Area law** (like CFT, unlike ETH) |
| **Quantifier** | Approximate KL conditions | Coherent information $I_c$ | Coherent information $I_c$ |
| **Key formula** | $\langle E_k|O|E_l\rangle \sim e^{-cN}$ | $\langle\phi_\beta|L_i|\phi_\alpha\rangle \sim n^{-\Delta}$ | $\langle S_m|\sigma|S_n\rangle \sim$ ? |

## Why Scars Might Form AQECC (Despite Violating ETH)

The Brandão et al. paper identifies two key conditions for AQECC:

1. **Diagonal condition:** $|\langle E_k|O|E_k\rangle - \langle E_{k'}|O|E_{k'}\rangle| \ll 1$ (local indistinguishability)
2. **Off-diagonal condition:** $|\langle E_k|O|E_l\rangle| \ll 1$ for $k \neq l$ (suppressed transitions)

**For ETH systems:** Both satisfied with exponential suppression $\sim e^{-cN}$.

**For Scar systems:** 
- Diagonal condition: Scar states have *similar* reduced density matrices (both area-law entangled)
- Off-diagonal condition: FSA algebra constrains transitions — only $Q^\pm$ connect scar levels

**Hypothesis:** The approximate $\mathfrak{su}(2)$ structure provides "selection rules" that suppress off-diagonal matrix elements for errors outside the algebra, analogous to how symmetry protects information.

## Code Space Definition

### Scar Encoding

Logical qubit encoded in Néel states:

$$|0_L\rangle = |\mathbb{Z}_2\rangle = |101010...\rangle$$
$$|1_L\rangle = |\mathbb{Z}_2'\rangle = |010101...\rangle$$

Reference-code entangled state for coherent information:

$$|\psi_{RQ}\rangle = \frac{1}{\sqrt{2}}\left(|0\rangle_R \otimes |\mathbb{Z}_2\rangle_Q + |1\rangle_R \otimes |\mathbb{Z}_2'\rangle_Q\right)$$

### Thermal Encoding (Comparison)

$$|\psi_{RQ}^{\rm th}\rangle = \frac{1}{\sqrt{2}}\left(|0\rangle_R \otimes |{\rm th}_1\rangle_Q + |1\rangle_R \otimes |{\rm th}_2\rangle_Q\right)$$

where $|{\rm th}_{1,2}\rangle$ are orthogonal thermal eigenstates at $E=0$.

## Noise Model

### Stage 1: Unitary Dynamics (Warm-up)

Channel: $\mathcal{N}_t(\rho) = e^{-iH_{\rm PXP}t} \rho \, e^{iH_{\rm PXP}t}$

- Scar encoding: $I_c(t)$ oscillates with revivals at $T_{\rm rev}$
- Thermal encoding: $I_c(t)$ decays monotonically

### Stage 2: Local Dephasing (Main Result)

Uniform dephasing channel:

$$\mathcal{N}_{p,\alpha}(\rho) = \bigotimes_{j=1}^{L} \left[(1-\frac{p}{2})\rho + \frac{p}{2}\sigma_\alpha^{[j]} \rho \, \sigma_\alpha^{[j]}\right]$$

for $\alpha \in \{X, Y, Z\}$.

### Key Observables

1. **Coherent information:** $I_c(R\rangle Q) = S(\rho_Q) - S(\rho_{RQ})$
2. **Threshold:** $p_c$ where $I_c \to 0$
3. **Scaling:** $p_c(L)$ behavior as $L \to \infty$

**Hypothesis:** $p_c^{\rm scar} > p_c^{\rm th}$ (or $p_c^{\rm scar} > 0$ while $p_c^{\rm th} = 0$)

## Numerical Protocol

### System Sizes
$L = 10, 12, 14, 16, 18, 20, 22$ (exact diagonalization)

### Computation Steps

1. **Prepare** $\rho_{RQ}^{(0)} = |\psi_{RQ}\rangle\langle\psi_{RQ}|$

2. **Apply noise:** $\rho_{RQ} = \mathcal{N}_{p,\alpha}(\rho_{RQ}^{(0)})$
   - Kraus operators for dephasing: $K_0 = \sqrt{1-p/2}\,I$, $K_1 = \sqrt{p/2}\,\sigma_\alpha$

3. **Compute reduced density matrices:**
   - $\rho_Q = {\rm Tr}_R(\rho_{RQ})$
   - $\rho_R = {\rm Tr}_Q(\rho_{RQ})$

4. **Compute entropies:** $S(\rho_{RQ})$, $S(\rho_Q)$

5. **Coherent information:** $I_c = S(\rho_Q) - S(\rho_{RQ})$

### Parameter Scans

- Noise strength: $p \in [0, 1]$, step $0.02$
- Noise type: $\alpha \in \{X, Y, Z\}$
- Stage 1 time: $t \in [0, 5T_{\rm rev}]$

## Theoretical Framework

### CFT ↔ ETH ↔ Scar Correspondence

| Concept | ETH Code | CFT Code | Scar Code |
|---------|----------|----------|-----------|
| Code subspace | $\{|E_k\rangle\}$ at energy density $\epsilon$ | $\{|\phi_\alpha\rangle\}$ low-energy | $\{|S_n\rangle\}$ scar tower |
| Algebraic structure | None (chaos) | Virasoro algebra | Approximate $\mathfrak{su}(2)$ (FSA) |
| KL diagonal | ETH: $\langle E_k|O|E_k\rangle \approx \bar{O}(\epsilon)$ | CFT: $\langle\phi|O|\phi\rangle \sim n^{-\Delta}$ | Scar: $\langle S_m|O|S_m\rangle \approx$ ? |
| KL off-diagonal | $|\langle E_k|O|E_l\rangle| \sim e^{-cN}$ | $|\langle\phi_\beta|O|\phi_\alpha\rangle| \sim n^{-\Delta}$ | $|\langle S_m|O|S_n\rangle| \sim$ ? |
| $b(n)$ coefficient | $\to 0$ (ETH) | $\propto n^{1-2\Delta}$ | TBD |
| Scaling exponent $\nu$ | $< 0$ always | $\nu = 1 - 2\Delta$ | TBD |
| Threshold exists? | Yes (any $p < 1$) | Iff $\Delta > 1/2$ | **To determine** |

### Key Analysis: Approximate Knill-Laflamme Conditions

Following Brandão et al., verify for error operators $E \in \{I, \sigma_\alpha^{[j]}\}$:

$$\langle \bar{i} | E | \bar{j} \rangle = C_E \delta_{ij} + \varepsilon_{ij}$$

**Compute for scar encoding:**
1. $\langle \mathbb{Z}_2 | \sigma_\alpha^{[j]} | \mathbb{Z}_2 \rangle$ — diagonal elements
2. $\langle \mathbb{Z}_2 | \sigma_\alpha^{[j]} | \mathbb{Z}_2' \rangle$ — off-diagonal elements  
3. $\varepsilon_{ij}$ scaling with system size $L$

**Key quantity (from Sang-Zou):**
$$b(L) = \frac{1}{D^2}\sum_{j=1}^{L}\left[D\,{\rm Tr}(\sigma_j P \sigma_j P) - {\rm Tr}(\sigma_j P){\rm Tr}(\sigma_j P)\right]$$

where $P = |\mathbb{Z}_2\rangle\langle\mathbb{Z}_2| + |\mathbb{Z}_2'\rangle\langle\mathbb{Z}_2'|$ is the code projector.

- If $b(L) \to 0$ as $L \to \infty$: code has threshold
- Scaling $b(L) \propto L^{1-2\Delta_{\rm eff}}$ determines effective "scar dimension"

## Expected Results

### Main Figures

1. **Coherent information vs noise strength**
   - $I_c(p)$ for scar vs thermal encoding
   - System sizes $L = 10, ..., 22$
   - Panels: $X$, $Y$, $Z$ dephasing

2. **Scaling collapse** (following Sang-Zou)
   - $I_c(p, L)$ vs $p \cdot L^\nu$ 
   - Extract effective exponent $\nu_{\rm scar}$
   - Compare with thermal states

3. **$b(L)$ coefficient**
   - Direct computation of KL violation
   - Scaling: $b(L) \propto L^{1-2\Delta_{\rm eff}}$?
   - Determine effective "scar scaling dimension"

4. **Threshold phase diagram**
   - $p_c(L)$ extrapolation to thermodynamic limit
   - Compare: $p_c^{\rm scar}$ vs $p_c^{\rm thermal}$ vs $p_c^{\rm CFT}$

5. **Dynamical coherent information**
   - $I_c(t)$ revivals for scar encoding
   - Decay for thermal encoding

### Key Claims

1. Scar states form AQECC despite violating ETH
2. Protection mechanism: approximate $\mathfrak{su}(2)$ algebra (distinct from ETH local indistinguishability)
3. Effective "scar scaling dimension" $\Delta_{\rm scar}$ determines threshold
4. Scars occupy a unique position: area-law entanglement (like CFT ground states) at finite energy density (like ETH states)

## Generalization

### Other Scar Systems

- AKLT chain (exact scars, SU(2) symmetry)
- XY model scars
- $\eta$-pairing in Hubbard model
- Rydberg ladders / 2D arrays

### Generalization Principle

Any system with approximate algebraic structure (tower of states with selection rules) can furnish a scar code. The "scaling dimension" analogue is the decay rate of matrix elements $\langle {\rm scar}_m | O_{\rm error} | {\rm scar}_n \rangle$ with system size.

### Experimental Connection

- Rydberg atom arrays: natural platform
- Realistic noise: dephasing, atom loss
- Coherent information measurable via tomography

## References

### Primary Reference
- **Brandão, Crosson, Şahinoğlu, Harrow**, "Quantum Error Correcting Codes in Eigenstates of Translation-Invariant Spin Chains," PRL 123, 110502 (2019), arXiv:1710.04631
  - Establishes ETH → AQECC connection
  - Provides approximate Knill-Laflamme framework
  - Shows translation-invariance alone yields AQECC

### Secondary Reference  
- **Sang & Zou**, "Approximate quantum error correcting codes from conformal field theory," arXiv:2406.09555
  - Coherent information methodology
  - Scaling collapse analysis
  - Threshold condition $\Delta > 1/2$

### Background
- Existing ergotropy paper (this project)
- FSA and scar tower literature (Papić et al.)
- Beny & Oreshkov, approximate QEC framework
