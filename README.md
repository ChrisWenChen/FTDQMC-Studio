# FTDQMC-Studio

Milestone 0: project skeleton, config parsing, logging, and placeholder pipeline stages.

## Build

Requirements:
- CMake >= 3.20
- C++20 compiler (`g++`, `clang++`, or `icpx`)
- MPI (default ON)
- Optional: PETSc / SLEPc / HDF5

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
```

MPI can be disabled for local smoke tests:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DFTDQMC_ENABLE_MPI=OFF
cmake --build build -j
```

## Run

CLI is fixed (no custom flags):

```
ftdqmc_studio <input_file> [PETSc/SLEPc options...]
```

Example (MPI):

```bash
mpirun -n 2 ./build/ftdqmc_studio docs/input_examples/minimal.toml
```

Example (PETSc/SLEPc options pass-through):

```bash
mpirun -n 2 ./build/ftdqmc_studio docs/input_examples/minimal.toml -eps_type krylovschur -eps_nev 4
```

Outputs are written to `[io].out_dir` from the TOML input. Per-rank logs are written by default, and rank 0 also writes an overview log.
