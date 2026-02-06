# FTDQMC-Studio — Project Specification (SPEC.md)

> **Single source of truth.** If any other document, code, or tool output conflicts with this file, **this file wins**.

## 0. Scope and Intent

FTDQMC-Studio is a general-purpose scientific computing software centered around Determinant Quantum Monte Carlo (DQMC) infrastructure and related linear-algebra workflows, designed for both local workstations and HPC systems.

**Primary near-term goal (pre-MC):**
- Generate/read auxiliary fields
- Build slice propagators \(B_\ell^\sigma\) and their chain products (B-chain)
- Build **time–space formulation** operator \(M_{mat}^\sigma\) (also referred to as “spacetime matrix”)
- Compute spectra: eigenvalues and singular values (dense or PETSc/SLEPc)

**Markov-chain Monte Carlo (MC) is explicitly NOT part of Milestone 0–5 unless stated in ROADMAP.md.** The project must not “accidentally” implement MC features early.

---

## 1. Hard Rules (Non-Negotiable)

### 1.1 Command Line Interface (CLI)
- The program accepts **no custom CLI flags** (no `--in`, no `--mode`, no `--help` beyond a simple usage message).
- Invocation is strictly:

  `ftdqmc_studio <input_file> [PETSc/SLEPc options ...]`

- `<input_file>` **must be `argv[1]`** and **must NOT start with `-`**.
  - If missing or invalid: print usage and exit with non-zero status.
- All arguments after the input file are reserved for PETSc/SLEPc and must be passed through without interpretation.

### 1.2 Input Configuration Format
- The only supported user configuration format is **TOML**.
- TOML is the authoritative source for:
  - workflow mode (spectrum vs mc)
  - lattice/model parameters
  - HS transform selection
  - field generation/reading
  - output settings
- PETSc/SLEPc solver options are configured **only** via PETSc/SLEPc command-line options.

### 1.3 Terminology and Definitions (Canonical)
- **B-slice propagator:** \(B_\ell^\sigma\) for imaginary-time slice \(\ell\), spin/flavor \(\sigma\).
- **B-chain:** the ordered product \( \mathcal{B}^\sigma = B_{Ltau-1}^\sigma B_{Ltau-2}^\sigma \cdots B_0^\sigma \).
- **Spacetime operator (time-space formulation):** \(M_{mat}^\sigma\), a block matrix of dimension \((N_s Ltau)\times(N_s Ltau)\) (with \(N_s\) spatial degrees of freedom, \(Ltau\) time slices).  
  - This is **NOT** the same as the conventional “\(M\)-matrix” that sometimes denotes \(I+\prod B\).  
- **Do not rename these terms** without updating SPEC.md.

### 1.4 Backends and Portability
- Must build on HPC with:
  - C++20
  - CMake ≥ 3.20
  - MPI
- Must support compilers: `icpx`, `g++`, `clang++`.
- Linear algebra backends must be selectable at build time:
  - MKL / Apple Accelerate / OpenBLAS (dense BLAS/LAPACK)
  - PETSc (sparse, vectors/matrices, KSP)
  - SLEPc (eigensolver/SVD)
- HDF5 output support must be optional.

**Hard constraint:** The codebase must not use Eigen. Dense matrices must be stored in a project-defined, column-major layout compatible with MKL/Accelerate/OpenBLAS.

---

## 2. Repository Structure (Canonical)

Top-level layout (minimum):
- `CMakeLists.txt`
- `cmake/` (optional, find modules/helpers)
- `include/ftdqmc/` (public headers)
- `src/` (library implementation)
- `apps/ftdqmc_studio/main.cpp` (main entry)
- `input_examples/` (TOML examples)
- `docs/` (design notes, math notes, roadmap, schema docs)
- `tests/` (unit and regression tests)

**Rule:** avoid putting substantial logic directly in `main.cpp`. Main should orchestrate: init → config → pipeline → finalize.

---

## 3. Build System Requirements

### 3.1 CMake Options (Expected)
Provide options (names may vary but must exist and be documented):
- `FTDQMC_ENABLE_MPI` (default ON)
- `FTDQMC_ENABLE_PETSC` (default OFF)
- `FTDQMC_ENABLE_SLEPC` (default OFF; requires PETSc)
- `FTDQMC_ENABLE_HDF5` (default OFF)
- `FTDQMC_LINALG_BACKEND` (string enum or equivalent):
  - `MKL`, `ACCELERATE`, `OPENBLAS`, `PETSC`

### 3.2 Dependency Handling
- Prefer `find_package(...)` when possible.
- If PETSc/SLEPc are enabled:
  - configure build using PETSc-provided CMake config or environment variables (`PETSC_DIR`, `PETSC_ARCH`) and document required setup.
- Third-party header-only libs (e.g., `toml++`, `spdlog`) may be vendored or fetched via CMake FetchContent.

---

## 4. Runtime Initialization Rules (MPI / PETSc / Logging)

### 4.1 MPI and PETSc Ownership
- When PETSc is enabled:
  - use `PetscInitialize`/`PetscFinalize` as the primary lifecycle.
  - do not double-initialize MPI.
- When PETSc is not enabled:
  - use `MPI_Init`/`MPI_Finalize`.

### 4.2 Output Directory Creation
- Output directory is defined by TOML: `[io].out_dir` (default `./out`).
- Only rank 0 creates directories; then `MPI_Barrier` ensures all ranks proceed safely.

### 4.3 Logging
- Must record (at minimum):
  - timestamp
  - MPI size, rank
  - build type, compiler ID/version
  - git commit hash (if available)
  - enabled backends (BLAS/LAPACK/PETSc/SLEPc/HDF5)
  - input file path
  - parsed workflow mode and task summary
- Logging policy:
  - rank 0 logs to stdout + file
  - per-rank log files are supported and enabled by default (`log_rank<r>.txt`)

---

## 5. Configuration Contract (TOML)

### 5.1 Workflow Selection
- TOML must specify:
  - `[workflow].mode = "spectrum" | "mc"`

### 5.2 Task Model
- Recommended schema supports a task list:
  - `[[tasks]]`
    - `kind = "eigs" | "svd" | "measure" | ...`
    - `target = "M_mat" | "Bprod" | ...` (must be explicit)
- In spectrum mode, tasks must be restricted to spectral tasks (eigs/svd).
- In mc mode, tasks may include measurement and related operations.

### 5.3 Lattice / Model / Imaginary Time
- Must define:
  - lattice type (at least: honeycomb)
  - boundary condition (OBC/PBC; twist may be added later)
  - model type (at least: Hubbard on honeycomb; Haldane-Hubbard)
  - imaginary time grid: `beta`, `dtau` (or `Ltau`)

### 5.4 HS Transform and Field I/O
- HS transform must be selectable by name:
  - `[hs].type = "<name>"`
- Field handling must support:
  - generate patterns: ferromagnetic / antiferromagnetic / random
  - read from file

---

## 6. Numerical Objects (Pre-MC)

### 6.1 Slice Propagators and B-Forms
- Must support B-slice construction forms (interface-level at minimum):
  - `KV`, `VK`, `KVK`
- B-chain refers strictly to the ordered product across time slices.

### 6.2 Spacetime Operator \(M_{mat}^\sigma\)
- Must support building \(M_{mat}^\sigma\) explicitly as:
  - dense matrix (for small systems)
  - PETSc Mat (for large systems)

---

## 7. Spectrum Computation Contract

### 7.1 Dense Mode
- For small systems, allow full eigenvalue decomposition / SVD via BLAS/LAPACK-based routines (no Eigen).
- Output eigenvalues/singular values in a stable, documented format.

### 7.2 PETSc/SLEPc Mode
- When enabled, spectral tasks must use SLEPc:
  - EPS for eigenvalues
  - SVD for singular values
- Solver configuration is controlled via PETSc/SLEPc command-line options, not TOML.

---

## 8. Output Contract

### 8.1 Text Output (Baseline)
At minimum, always produce:
- `out_dir/run_summary.txt` (or equivalent)
  - includes config summary and success status

### 8.2 HDF5 Output (Optional)
If HDF5 enabled, use a structured layout (exact paths can evolve but must remain documented):
- `/meta` (build/runtime info)
- `/params` (config snapshot)
- `/spectrum/<task_name>/...` (eigs/svd results)
- `/field/...` (optional snapshots)

---

## 9. Testing Contract (Minimal)

### 9.1 Test Framework
- Use GoogleTest (gtest) for unit and regression tests.
- All tests must be written using gtest macros.

### 9.1 Build Tests
- CI or local tests must compile with:
  - g++
  - clang++
  - optional icpx (if available)

### 9.2 Smoke Tests
- A minimal example TOML must run end-to-end under MPI (even if algorithm steps are placeholders in Milestone 0).

### 9.3 Math Consistency (Later Milestones)
When \(\hat O^\sigma\) is implemented, add a small-system check:
- `det(M_mat)` must match `det(I + B_chain)` in dense mode (within tolerance), for small sizes.

---

## 10. Development Discipline

- Keep interfaces stable once introduced; deprecate explicitly if changes are necessary.
- Avoid “quick hacks” in `main.cpp`.
- Do not add custom CLI options.
- Do not introduce alternative config formats (no JSON, no YAML, no ad-hoc `.in`).
- Always update docs (`SPEC.md`, `ROADMAP.md`, `input_examples/*`) when behavior changes.

---

## 11. Versioning and Reproducibility

- Record git hash and build flags in outputs.
- Favor deterministic RNG seeding rules (when MC is introduced).
- Output must include sufficient metadata to reproduce results (backend selection, solver options summary if possible).

---

**End of SPEC.md**
