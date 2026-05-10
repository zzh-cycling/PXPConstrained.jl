# Report: Why `iso_K2MSS` fails in the `k = π` sectors

## Scope

I investigated the `k = π` (`k = L/2`) inversion-sector construction in:

- `src/PXPSymmetry.jl`:
  - `iso_K2MSS`
  - downstream `PXP_MSS_Ham` sector splitting logic

The issue was reproduced using direct K-space inversion diagnostics.

---

## Symptom

You observed:

- `PXP_MSS_Ham(12, 6, -1)` is only `5×5`
- expected scar states in `k=π, I=-1` are missing

This is consistent with an incorrect sector decomposition.

---

## Key finding

At `k=π`, inversion in the K basis is **not** just a plain swap `n ↔ nR` with unit sign.

For a representative `n`, if reflection maps to `nR` with translation offset `dR`,

\[
\mathcal I \, |n\rangle_k = \eta_n \, |n_R\rangle_k,\quad \eta_n = \omega_k^{d_R}.
\]

At `k=π`, \(\omega_k=-1\), so \(\eta_n = (-1)^{d_R}\), i.e. **state-dependent ±1 phase**.

Current `iso_K2MSS` ignores this phase and assumes fixed `±` combinations.  
It also drops all `n==nR` states in `inv=-1`, but for `k=π` many such states have \(\eta_n=-1\) and must belong to the `I=-1` sector.

---

## Evidence

I built the K-space inversion operator \(I_k\) directly and compared:

1. **True sector dimensions** from eigenspaces of \(I_k\)
2. **Current `iso_K2MSS` dimensions**
3. Inversion eigenvalue content of current `iso_K2MSS` columns (`W' I_k W`)

### Results

| L | K dim | self \(\eta=+1\) | self \(\eta=-1\) | non-self pairs | expected dims (+,-) | true eig dims (+,-) | current iso dims (+,-) |
|---|---:|---:|---:|---:|---:|---:|---:|
| 12 | 29 | 3 | 16 | 5 | (8, 21) | (8, 21) | (24, 5) |
| 14 | 59 | 4 | 25 | 15 | (19, 40) | (19, 40) | (44, 15) |
| 16 | 142 | 10 | 44 | 44 | (54, 88) | (54, 88) | (98, 44) |

And current `iso_K2MSS` subspaces are mixed:

- For `L=12`, `iso_K2MSS(..., inv=1)` columns contain 6 true `+1` and 18 true `-1` eigenvectors.
- For `L=12`, `iso_K2MSS(..., inv=-1)` columns contain 2 true `+1` and 3 true `-1`.

So the current `inv` labels do **not** match true inversion sectors at `k=π`.

---

## Root cause (precise)

In `iso_K2MSS`, the basis construction assumes:

- `inv=+1` ⇒ symmetric combination of paired representatives
- `inv=-1` ⇒ antisymmetric combination
- `n==nR` always belongs to `inv=+1`

This is valid for `k=0` (phase is always +1), but invalid for `k=π` where \(\eta_n=(-1)^{d_R}\) must be included.

---

## Correct first-principles sector rule at `k=π`

Let \(\eta_n = (-1)^{d_R(n)}\).

1. **Self-reflection case** (`n==nR`):
   - state belongs to inversion eigenvalue \(\lambda=\eta_n\) (can be +1 or -1)

2. **Paired case** (`n \neq nR`):
   - inversion eigenstates are
   \[
   |\psi_\lambda\rangle \propto |n\rangle_k + \lambda\,\eta_n\,|n_R\rangle_k,\quad \lambda=\pm1.
   \]

So each non-self pair contributes one state to each sector, and self states split by \(\eta_n\).

---

## Impact

- `iso_K2MSS(L, L/2, ±1)` currently defines mislabeled/mixed subspaces.
- Any `PXP_MSS_Ham` built on this split inherits incorrect sector dimensions and spectra.
- Missing expected scar states in the nominal `I=-1` block is a direct consequence.

---

## Recommendation

Revise `iso_K2MSS` (and consistent helpers: `PXP_MSS_basis`, `mapstate_MSS2K`, `PXP_MSS_Ham`) to include the inversion phase \(\eta_n\) at `k=π` exactly as above.

This is the minimal physics-correct fix and should restore correct `k=π` sector dimensions and scar visibility.
