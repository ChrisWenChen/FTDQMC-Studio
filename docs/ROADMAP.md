# ROADMAP.md

## Goal

This roadmap fixes the **milestone order** for the current development phase:
**field generation/IO → model/lattice → build slice propagators → construct time-space `M_mat^σ` → compute spectra**.

Hard rule:
- **Do not implement any Markov-chain / Monte Carlo updates** (no local flips, no accept/reject, no measurements requiring sampling) until explicitly scheduled in a later roadmap revision.

---

## Milestone 0 — Project Skeleton (Build / Config / Logging)

Deliverables:
- C++20 + CMake project that builds with `g++`, `clang++`, optional `icpx`.
- MPI runnable (`mpirun -n ...`).
- CLI convention (no custom flags):
  - `ftdqmc_studio <input_file> [PETSc/SLEPc options...]`
  - only `<input_file>` is handled by our program; all other args are for PETSc/SLEPc.
- TOML config parsing (minimal fields).
- Logging with rank-aware log files.
- Output directory creation (rank0 + barrier).
- Placeholders for pipeline stages:
  - “field: placeholder”
  - “build M_mat^σ: placeholder”
  - “spectrum: placeholder”
  - “write output: dummy file”

---

## Milestone 1 — Lattice + K (Honeycomb Hubbard + Haldane-Hubbard)

Deliverables:
- Honeycomb lattice construction:
  - geometry, site indexing (time-major linearization), sublattice labeling, neighbor lists.
  - boundary conditions (at least OBC/PBC; twist BC can be an interface stub).
- Model definitions:
  - Hubbard onsite parameters (U, μ).
  - Haldane terms: NN `t1`, NNN `t2`, complex phase `φ`, and optional staggered mass.
- Kinetic matrix builder `K`:
  - supports complex hoppings (required for Haldane-Hubbard).
  - backend-agnostic interface (dense now, PETSc later).

---

## Milestone 2 — HS Plugin + Field Generation/Read/Write

Deliverables:
- HS transformation as a plugin interface:
  - selectable HS type via input TOML.
  - provides mapping from auxiliary field `s(i,l)` to diagonal potentials `V_l^σ`.
  - **Clarification:** $V$ does not include the $\Delta\tau$ factor (applied by propagator builder).
- Auxiliary field workflows (no sampling):
  - generate: ferromagnetic / antiferromagnetic / random patterns
  - read from file
  - write to file
- IO formats:
  - at least text
  - HDF5 optional (if enabled)

---

## Milestone 3 — Slice Propagators and B-chain (KV / VK / KVK)

Deliverables:
- Build `B_l^σ` for each slice and spin using selected B-form:
  - KV, VK, KVK (symmetric).
- Define and compute/apply the ordered product (0 to $Ltau-1$):
  - `Bchain^σ = B_{Ltau-1}^σ ... B_0^σ`
- Provide both:
  - explicit matrix form (dense baseline)
  - apply-to-vector form (for future matrix-free mode)

No Markov-chain updates.

---

## Milestone 4 — Explicit Construction of Time-space `M_mat^σ` (dense + PETSc)

Deliverables:
- Construct time-space block matrix `M_mat^σ` with dimension `(Ns*Ltau) × (Ns*Ltau)`,
  matching the fixed block/sign conventions (lower-subdiagonal with wrap-around) in `MATH_NOTES.md`.
- Backends:
  - dense explicit (baseline)
  - PETSc `Mat` explicit (optional, if enabled)
  - **MPI Strategy:** Row-wise distribution of blocks; each rank owns a subset of time slices.

Must include regression tests:
- block placement/sign test on small cases
- determinant identity:
  - `det(M_mat^σ) == det(I + Bchain^σ)` for small sizes (dense baseline)
  - **Stability:** Use QR-based stabilization for the `det(I + Bchain)` side.

---

## Milestone 5 — Spectra of `M_mat^σ` (Eigenvalues and Singular Values)

Deliverables:
- Compute spectra per spin/flavor:
  - eigenvalues (full diagonalization or iterative)
  - singular values (full SVD or iterative)
- Backend options:
  - dense (LAPACK/Eigen) for full methods
  - SLEPc (EPS/SVD) for iterative methods and PETSc matrices
- **Diagnostics:** Optional calculation of the condition number of `M_mat^σ` (testing only).
- Output:
  - write eigenvalues/singular values + metadata to text and/or HDF5

---

## Future (Not in scope yet): Markov-chain Monte Carlo

Only after Milestone 5 is complete and validated:
- local updates, acceptance ratios, wrapping, stabilization stacks, measurements, binning, etc.
