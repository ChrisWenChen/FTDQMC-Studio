# FTDQMC-Studio — Architecture (ARCHITECTURE.md)

This document defines the canonical repository layout, module responsibilities, dependency directions, and the stable public interfaces.  
If any code generation diverges from this architecture, refactor back to match this file.

---

## 0. Architecture Goals

- **HPC-first**: MPI-ready, PETSc/SLEPc-friendly, scalable I/O.
- **No CLI drift**: configuration comes from a single TOML file; PETSc/SLEPc options remain in argv untouched (see SPEC.md).
- **Clear layering**: separate *configuration*, *domain model (lattice/model/HS)*, *numerical kernels*, and *runners (workflows)*.
- **Backend isolation**: dense (BLAS/LAPACK/Eigen if allowed) and PETSc/SLEPc implementations are swappable behind stable interfaces.
- **Algorithm staging**: early milestones focus on field generation/reading, spacetime operator construction, and spectrum tasks; MC is added later without rewriting core modules.
- **Parallelism**: In MPI-enabled PETSc mode, the `M_mat` is distributed such that each rank is responsible for a subset of time slices (contiguous row blocks).

---

## 1. High-Level Layers

### 1.1 Layers and Responsibilities

1) **App Layer**
- `apps/ftdqmc_studio/main.cpp`
- Owns: initialization, configuration loading, selecting workflow runner, top-level error handling.
- Must contain minimal logic: no math kernels.

2) **Core Layer**
- `src/core`, `include/ftdqmc/core`
- Owns: MPI/PETSc lifecycle wrapper, logging, filesystem/output directory, error utilities, small helpers.

3) **Config Layer**
- `src/config`, `include/ftdqmc/config`
- Owns: TOML parsing and validation into strong typed structs.
- Converts string enums into `enum class` values early.

4) **Domain Layer (Physics/Modeling)**
- `src/lattice`, `src/model`, `src/hs`, `src/field`, `include/ftdqmc/...`
- Owns: lattice geometry, boundary conditions, model parameters, HS transforms, auxiliary field objects.

5) **Numerics Layer**
- `src/linalg`, `src/propagator`, `src/spacetime`, `src/spectrum`, `include/ftdqmc/...`
- Owns: dense/PETSc matrix representations, slice propagator builders (KV/VK/KVK), B-chain application, spacetime operator (`M_mat`) construction, eigen/SVD tasks.

6) **Workflow Layer**
- `src/workflow`, `include/ftdqmc/workflow`
- Owns: `SpectrumRunner` and later `MCRunner`.
- Orchestrates pipelines using Domain + Numerics modules.

7) **I/O Layer**
- `src/io`, `include/ftdqmc/io`
- Owns: text and HDF5 output format, structured metadata writing, snapshots (field), and results writing (spectrum).

---

## 2. Dependency Rules (Directed Graph)

### 2.1 Allowed Dependencies
- `apps` → `workflow`, `core`, `config`
- `workflow` → `domain`, `numerics`, `io`, `core`
- `numerics` → `domain` (read-only), `linalg`, `core`
- `domain` → `core` only (for logging/error helpers, no numerics)
- `io` → `core`, `config` (for metadata), and minimal types from `domain`/`numerics` as needed

### 2.2 Forbidden Dependencies
- `domain` must not depend on `numerics` or `workflow`.
- `config` must not depend on `numerics` kernels (parsing must stay lightweight).
- `io` must not call solvers or builders; it only serializes results and metadata.
- `main.cpp` must not contain substantial numerical logic.

---

## 3. Canonical Directory Layout

├── CMakeLists.txt
├── cmake/
│ ├── Modules/ # optional: FindMKL.cmake, etc.
│ └── Toolchains/ # optional: HPC toolchain presets
├── include/
│ └── ftdqmc/
│ ├── core/
│ ├── config/
│ ├── lattice/
│ ├── model/
│ ├── hs/
│ ├── field/
│ ├── linalg/
│ ├── propagator/
│ ├── spacetime/
│ ├── spectrum/
│ ├── io/
│ └── workflow/
├── src/
│ ├── core/
│ ├── config/
│ ├── lattice/
│ ├── model/
│ ├── hs/
│ ├── field/
│ ├── linalg/
│ ├── propagator/
│ ├── spacetime/
│ ├── spectrum/
│ ├── io/
│ └── workflow/
├── apps/
│ └── ftdqmc_studio/
│ └── main.cpp
├── input_examples/
│ ├── minimal.toml
│ ├── spectrum_honeycomb.toml
│ └── spectrum_haldane_hubbard.toml
├── docs/
│ ├── SPEC.md
│ ├── ARCHITECTURE.md
│ ├── ROADMAP.md
│ └── MATH_NOTES.md
└── tests/
├── unit/
└── regression/


---

## 4. Canonical Names and Key Concepts

### 4.1 Core Symbols (Must Match SPEC.md)
- `B_l^sigma` is the slice propagator.
- `B-chain` is `B_{Ltau-1} ... B_0` (ordered product across time slices).
- `M_mat` is the **spacetime/time–space** block operator used for spectral computations.
- Avoid ambiguous naming like “M matrix” without qualifier. In code, prefer:
  - `SpacetimeOperator` / `M_mat`
  - `BChain` for the product and/or apply operation
  - `Mprod = I + B_chain` if later needed

### 4.2 Indexing Convention (Canonical)
- Time slice index: `l = 0..Ltau-1`
- Spatial site index: `i = 0..Ns-1`
- Spin/flavor index: `sigma` enumerated (e.g., Up/Down)
- B-chain order: `B_{Ltau-1} ... B_0` (rightmost acts first on vectors)
- Linearization: `idx = l * Ns + i` (time-major)

---

## 5. Public API Contracts (Headers)

This section defines stable interface “shapes”. Actual implementation details may evolve, but these class/function responsibilities must remain.

### 5.1 Core: Environment, Logging, Error

**Header:** `include/ftdqmc/core/environment.hpp`
- `struct Environment`
  - owns MPI communicator info (rank/size) for non-PETSc mode
  - for PETSc mode, provides wrappers to query rank/size. Each rank handles a subset of time slices.
  - stores output directory (resolved path)
  - stores backend flags and build info (best-effort)

**Header:** `include/ftdqmc/core/logging.hpp`
- `void init_logging(const Environment&, const std::string& out_dir);`
- `Logger& log();` (or equivalent)

**Header:** `include/ftdqmc/core/status.hpp`
- `struct Status` for run summary and error codes

### 5.2 Config: TOML → Strong Types

**Header:** `include/ftdqmc/config/config.hpp`
- `enum class WorkflowMode { Spectrum, MC };`
- `enum class TaskKind { Eigs, Svd, Measure };`
- `enum class TargetKind { M_mat, Bprod, Bchain };` (expandable)
- `struct Config`
  - `WorkflowMode workflow_mode;`
  - lattice/model/imag_time/hs/field/io sections
  - `std::vector<TaskSpec> tasks;`
- `Config load_config_from_toml(const std::string& path);`
- `void validate_config(const Config&);` (throws or returns Status)

**Rule:** Convert all string enums to `enum class` during load.

### 5.3 Lattice

**Header:** `include/ftdqmc/lattice/lattice.hpp`
- `enum class LatticeType { Honeycomb, Chain1D, Rectangular, Hexagonal };`
- `enum class BoundaryCondition { OBC, PBC, TBC };`
- `struct LatticeSpec { LatticeType type; int Lx, Ly; BoundaryCondition bc; /* twist params */ };`
- `class Lattice`
  - `int num_sites() const;`
  - `std::vector<int> neighbors(int site) const;` (or adjacency list)
  - `/* optional */ sublattice(site), coordinates(site)`

**Honeycomb specialization**
- `Lattice make_honeycomb(const LatticeSpec&);`

### 5.4 Model and K Construction

**Header:** `include/ftdqmc/model/model.hpp`
- `enum class ModelType { Hubbard, HaldaneHubbard };`
- `struct ModelSpec { ModelType type; double U; double mu; /* Haldane params: t1,t2,phi,staggered... */ };`

**Header:** `include/ftdqmc/model/kbuilder.hpp`
- `class KBuilder`
  - `DenseOrPetscMatrix build_K(const Lattice&, const ModelSpec&, const ImagTimeSpec&);`
  - `DenseOrPetscMatrix build_expK(const DenseOrPetscMatrix& K, double dtau);` (optional; may move to numerics)

### 5.5 Imaginary Time

**Header:** `include/ftdqmc/core/imag_time.hpp`
- `struct ImagTimeSpec { double beta; double dtau; int Ltau; };`
- Canonical relation: `Ltau = round(beta / dtau)` unless specified explicitly.

### 5.6 HS Transform (Pluggable)

**Header:** `include/ftdqmc/hs/hs.hpp`
- `enum class HSType { HirschDiscrete, /* future */ };`
- `struct HSSpec { HSType type; /* channel, lambda policy, etc */ };`

**Header:** `include/ftdqmc/hs/transform.hpp`
- `class IHSTransform`
  - virtual destructor
  - `int num_flavors() const;`
  - `FieldAlphabet alphabet() const;`  (e.g., Ising ±1)
  - `void build_diag_potential(/*out*/ std::vector<double or complex>& v,
                              int sigma, int l,
                              gsl::span<const int or double> field_slice,
                              const HSSpec&, const ModelSpec&, const ImagTimeSpec&) const;`
- `std::unique_ptr<IHSTransform> make_hs_transform(const HSSpec&);`

**Rule:** HS must not depend on solver backends; it produces diagonal potentials only.

### 5.7 Auxiliary Field

**Header:** `include/ftdqmc/field/field.hpp`
- `enum class FieldInit { Ferromagnetic, Antiferromagnetic, Random, FromFile };`
- `struct FieldSpec { FieldInit init; std::string file; uint64_t seed; };`
- `class AuxiliaryField`
  - stores `field[i, l]` (layout must be documented; e.g., contiguous by time)
  - `int ns() const; int ltau() const;`
  - `span<const int> slice(int l) const;`
  - `span<int> slice(int l);`

**Header:** `include/ftdqmc/field/io.hpp`
- `AuxiliaryField read_field(const std::string& path);`
- `void write_field(const std::string& path, const AuxiliaryField&);`
- Must support text and HDF5 (if enabled).

### 5.8 Linear Algebra Abstractions

**Header:** `include/ftdqmc/linalg/matrix.hpp`
- Provide minimal unified types:
  - `DenseMatrix` (complex or real)
  - `PetscMatWrapper` (when PETSc enabled)
- Provide a tagged union or variant:
  - `using Matrix = std::variant<DenseMatrix, PetscMatWrapper>;`

**Rule:** Keep the “variant boundary” narrow; most kernels should be specialized by backend.

### 5.9 Slice Propagators and B-Chain

**Header:** `include/ftdqmc/propagator/builder.hpp`
- `enum class BForm { KV, VK, KVK };`
- `struct PropagatorSpec { BForm form; };`
- `class BSliceBuilder`
  - `Matrix build_B_slice(int sigma, int l,
                          const Matrix& expK,
                          const AuxiliaryField& field,
                          const IHSTransform& hs,
                          const PropagatorSpec&,
                          const Config& /* or needed specs */);`
  - `Matrix build_expV(/* diag potential */);` (may be internal)

**Header:** `include/ftdqmc/propagator/bchain.hpp`
- `class BChain`
  - stores per-slice `B_l^sigma` or provides apply-only access
  - `Matrix product(int sigma) const;` (dense mode)
  - `void apply(int sigma, /*in*/ Vec x, /*out*/ Vec y) const;` (operator mode; optional later)

### 5.10 Spacetime Operator (M_mat)

**Header:** `include/ftdqmc/spacetime/m_mat.hpp`
- `class SpacetimeOperatorBuilder`
  - `Matrix build_M_mat(int sigma, const std::vector<Matrix>& B_slices,
                       int ns, int ltau);`
- Must support:
  - dense explicit construction
  - PETSc Mat explicit block construction

### 5.11 Spectrum Tasks

**Header:** `include/ftdqmc/spectrum/spectrum.hpp`
- `enum class SpectrumKind { Eigs, Svd };`
- `enum class SolverBackend { Dense, SLEPc };`
- `struct SpectrumTaskSpec`
  - kind, target, solver backend, nev, tol, max_it, which, etc
- `struct SpectrumResult`
  - eigenvalues or singular values (vector of real/complex as appropriate)
  - optional residuals / iteration stats
- `class SpectrumSolver`
  - `SpectrumResult run(const SpectrumTaskSpec&, const Matrix& target_matrix, const Environment&);`

**Rule:** For SLEPc mode, all solver fine-tuning comes from PETSc options (argv).

### 5.12 Workflow Runners

**Header:** `include/ftdqmc/workflow/runner.hpp`
- `class IRunner { virtual Status run() = 0; };`
- `std::unique_ptr<IRunner> make_runner(const Config&, const Environment&);`

**Header:** `include/ftdqmc/workflow/spectrum_runner.hpp`
- `class SpectrumRunner : public IRunner`
  - pipeline:
    1) lattice/model/time specs
    2) HS transform
    3) field generate/read
    4) build expK
    5) build B slices
    6) build M_mat
    7) execute tasks (eigs/svd)
    8) write outputs

**Header:** `include/ftdqmc/workflow/mc_runner.hpp` (placeholder until MC milestone)
- defined but not implemented initially; must not block spectrum work.

---

## 6. Data Flow (Spectrum Mode)

In `workflow.mode = "spectrum"`:

1) `Config` loaded (TOML) → validated
2) `Environment` initialized (MPI/PETSc) → logging ready → out_dir created
3) Domain setup:
   - `Lattice` built from `LatticeSpec`
   - `ModelSpec` chosen (Hubbard / Haldane-Hubbard)
   - `ImagTimeSpec` computed
   - `IHSTransform` created
   - `AuxiliaryField` generated or read
4) Numerics setup:
   - `K` built (dense or PETSc)
   - `expK = exp(-dtau * K)` (dense or PETSc strategy)
   - B slices built for each `l` and `sigma`
   - `M_mat` built for each `sigma`
5) Spectrum tasks executed:
   - each `[[tasks]]` consumes the chosen `target`
6) Output:
   - always write a summary text file
   - write spectrum arrays (text or HDF5 if enabled)
   - optionally write field snapshot and metadata

---

## 7. Output Organization (Canonical)

Inside `[io].out_dir`:

- `run_summary.txt` (always)
- `logs/`
  - `log_rank0.txt`
  - `log_rank<r>.txt` (optional but default on)
- `spectrum/`
  - `task_<idx>_<kind>_<target>_sigma<up|down>.txt` (baseline)
  - `task_<idx>_<kind>_<target>.h5` (if HDF5 enabled)
- `field/`
  - `field_snapshot.txt` or `field_snapshot.h5` (optional)

---

## 8. Extension Points (Planned Growth)

- Add more lattices: chain1D, rectangular, hexagonal (graph + metric)
- Add more HS transforms (SU(2) symmetric, etc.)
- Add operator-mode M_mat/B-chain application (matrix-free)
- Add stabilization stacks for B-chain (QR/UDT/SVD) when transitioning to MC
- Add MC workflow without changing:
  - config contract
  - CLI contract
  - module boundaries

---

## 9. “Do Not Drift” Checklist

When generating or refactoring code, verify:

- No custom CLI flags were introduced.
- TOML remains the only config file.
- `M_mat` is consistently used for spacetime operator naming.
- Domain modules do not depend on Numerics/Workflow.
- `main.cpp` stays thin.
- PETSc/SLEPc options pass-through remains intact.

---

**End of ARCHITECTURE.md**
