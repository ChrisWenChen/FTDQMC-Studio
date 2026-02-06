# MATH_NOTES.md

## Purpose

This note fixes the **mathematical definitions and naming conventions** for the current development target:
**time-space (space–imaginary-time) formulation** and its associated **block matrix `M_mat^σ`**.

Hard rule for the codebase:
- **`M_mat`** (or `M_mat_sigma`) refers **ONLY** to the **time-space block matrix** described below.
- The conventional equal-time object `I + Π B` (sometimes also called "M" in other literature) must be named explicitly as:
  - `I_plus_Bchain` or `M_equal_time` (choose one and stick to it),
  - but **never** call it `M_mat`.

This prevents symbol drift and wrong sign/block placements during implementation.

---

## Core Objects and Dimensions

Let:

- `Ns` = number of spatial lattice sites (single-particle basis size per spin).
- `Ltau` = number of imaginary-time slices (total index $l = 0, 1, \dots, Ltau-1$).
- `σ`  = spin/flavor index (e.g., `up`, `down`).
- `B_l^σ` = slice propagator (an `Ns × Ns` matrix) for time slice `l = 0..Ltau-1`.

### Linearization Formula

For the time-space block matrix, a global index $I$ is mapped from space index $i \in [0, Ns-1]$ and time index $l \in [0, Ltau-1]$ as:
- **$I = l \cdot Ns + i$** (Time-major ordering)

### Time-space matrix `M_mat^σ`

Define `M_mat^σ` as an `(Ns*Ltau) × (Ns*Ltau)` block matrix with `Ltau × Ltau` blocks, each block of size `Ns × Ns`.

So the total dimension is:

- **`M_mat^σ ∈ C^{(Ns*Ltau) × (Ns*Ltau)}`**

(Real/complex depends on the model; Haldane-Hubbard generally leads to complex hoppings.)

---

## Block Structure of `M_mat^σ`

### Block layout (Lower-triangular cyclic)

`M_mat^σ` is block-cyclic with:

- Identity blocks on the diagonal,
- `-B_l^σ` on the first sub-diagonal,
- and a wrap-around block `+B_{Ltau-1}^σ` in the top-right corner.

In block form (0-indexed slices):

M_mat^σ =
| I            0            0      ...  +B_{Ltau-1}^σ |
| -B_0^σ       I            0      ...  0             |
| 0            -B_1^σ       I      ...  0             |
| :            :            :      ...  :             |
| 0            0            0      ...  I             |


### Indexing convention

We use **0-based indexing** consistently:

- time slice blocks are ordered `0, 1, ..., Ltau-1`
- diagonal blocks: $(l, l) \implies I$ (Identity)
- sub-diagonal blocks: $(l+1, l) \implies -B_l^σ$ for $l = 0 \dots Ltau-2$
- wrap-around block: $(0, Ltau-1) \implies +B_{Ltau-1}^σ$

Keep these signs and placements exactly. Do not “simplify” them.

---

## Relationship to the B-chain product

Define the **B-chain** (time-ordered product):

- **`Bchain^σ = B_{Ltau-1}^σ · B_{Ltau-2}^σ · ... · B_0^σ`**

Then `M_mat^σ` satisfies the determinant identity:

- **`det(M_mat^σ) = det(I + Bchain^σ)`**

### Stability and QR Decomposition

For large $Ltau \cdot \beta$, the direct product `Bchain` and its determinant $\det(I + Bchain)$ are numerically unstable.
Validation tests should use a **QR-based stable determinant** calculation (similar to the UdV stabilization in DQMC) even for the equal-time identity check.

### Numerical Diagnostics

For testing and debugging, the **Condition Number** of `M_mat^σ` can be computed to assess the stability of the time–space formulation. This should be an optional task and disabled in production runs.

---

## Definition of Slice Propagators `B_l^σ`

The slice propagator is assembled from a kinetic part `K` and a potential/HS part `V_l^σ`.

We denote:
- `Δτ` = imaginary-time step.
- `K` = kinetic (hopping) matrix, `Ns × Ns`.
- `V_l^σ` = diagonal matrix derived from HS field at slice `l` for spin `σ`, `Ns × Ns`.

**Note:** $V_l^\sigma$ is the potential contribution and does **not** include the $\Delta\tau$ factor. The factor is applied during propagator construction.

---

## Supported B-forms (KV / VK / KVK)

We support the following constructions for each slice `l` and spin `σ`.

### 1) KV form

**Definition:**
- `B_l^σ = exp(-Δτ K) · exp(-Δτ V_l^σ)`

---

### 2) VK form

**Definition:**
- `B_l^σ = exp(-Δτ V_l^σ) · exp(-Δτ K)`

---

### 3) KVK (symmetric/Trotter-symmetric) form

**Definition:**
- `B_l^σ = exp(-(Δτ/2) K) · exp(-Δτ V_l^σ) · exp(-(Δτ/2) K)`

---

## HS Field and `V_l^σ`

We treat the HS transformation as a plugin that maps the auxiliary field `s(i,l)` to diagonal potentials for each spin/flavor:

- `V_l^σ = diag( v_1^σ(l), v_2^σ(l), ..., v_{Ns}^σ(l) )`

The HS plugin must provide:
- `v_i^σ(l)` for each site and slice. The builder then computes `exp(-Δτ v_i^σ(l))`.

---

## Implementation Guidance (to avoid symbol drift)

1) **Never call `I + Bchain` “M_mat”.**
   - `M_mat` is time-space block matrix.
   - `I + Bchain` is equal-time object.

2) **Parallelism Strategy:**
   - In MPI-enabled PETSc mode, the `M_mat` is distributed such that **each rank is responsible for a subset of time slices** (contiguous rows of blocks).

3) Validation tests must include:
   - dimension checks: `(Ns*Ltau) × (Ns*Ltau)`
   - block placement checks for a small `(Ns, Ltau)` case
   - determinant identity: `det(M_mat^σ) == det(I + Bchain^σ)` using QR-based stability for the right-hand side.

---

## Minimal Glossary

- **Ns**: number of spatial sites.
- **Ltau**: number of imaginary-time slices.
- **B_l^σ**: slice propagator for slice `l`.
- **Bchain^σ**: ordered product `B_{Ltau-1}^σ ... B_0^σ`.
- **M_mat^σ**: time-space block matrix of dimension `(Ns*Ltau) × (Ns*Ltau)`.
- **I_plus_Bchain**: `I + Bchain^σ` (equal-time object).
