#include <iostream>
#include <vector>
#include <cmath>

#ifdef FTDQMC_USE_OPENMP
  #include <omp.h>
#endif

#ifdef FTDQMC_USE_MPI
  #include <mpi.h>
#endif

#ifdef FTDQMC_USE_PETSC
  #include <petsc.h>
#endif

#ifdef FTDQMC_USE_SLEPC
  #include <slepceps.h>
#endif

// BLAS: Accelerate on macOS, cblas elsewhere
#ifdef __APPLE__
  #include <Accelerate/Accelerate.h>
#else
  #include <cblas.h>
#endif

static void print_build_flags() {
  std::cout << "FTDQMC link smoke test\n";

#if defined(FTDQMC_BLAS_MKL)
  std::cout << "  BLAS backend: MKL\n";
#elif defined(FTDQMC_BLAS_OPENBLAS)
  std::cout << "  BLAS backend: OpenBLAS\n";
#elif defined(FTDQMC_BLAS_ACCELERATE)
  std::cout << "  BLAS backend: Accelerate\n";
#else
  std::cout << "  BLAS backend: (unknown)\n";
#endif

#ifdef FTDQMC_USE_MPI
  std::cout << "  MPI: ON\n";
#else
  std::cout << "  MPI: OFF\n";
#endif

#ifdef FTDQMC_USE_OPENMP
  std::cout << "  OpenMP: ON\n";
#else
  std::cout << "  OpenMP: OFF\n";
#endif

#ifdef FTDQMC_USE_PETSC
  std::cout << "  PETSc: ON\n";
#else
  std::cout << "  PETSc: OFF\n";
#endif

#ifdef FTDQMC_USE_SLEPC
  std::cout << "  SLEPc: ON\n";
#else
  std::cout << "  SLEPc: OFF\n";
#endif
}

static void test_blas_dgemm() {
  const int n = 64;
  std::vector<double> A(n*n), B(n*n), C(n*n, 0.0);

  for (int i = 0; i < n*n; ++i) {
    A[i] = std::sin(0.001 * i);
    B[i] = std::cos(0.001 * i);
  }

  // C = A * B
  cblas_dgemm(
    CblasColMajor, CblasNoTrans, CblasNoTrans,
    n, n, n,
    1.0, A.data(), n,
         B.data(), n,
    0.0, C.data(), n
  );

  double checksum = 0.0;
  for (double v : C) checksum += v;

  std::cout << "  BLAS dgemm checksum = " << checksum << "\n";
}

#ifdef FTDQMC_USE_MPI
static void test_mpi_allreduce() {
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  double local = 1.0 + rank;
  double global = 0.0;
  MPI_Allreduce(&local, &global, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  if (rank == 0) {
    std::cout << "  MPI world size = " << size << "\n";
    std::cout << "  MPI Allreduce sum(1+rank) = " << global << "\n";
  }
}
#endif


#ifdef FTDQMC_USE_PETSC
static void test_petsc_vec_mat() {
  // Minimal PETSc objects to ensure link + basic functionality
  Vec x, y;
  Mat A;

  const PetscInt n = 10;

  VecCreate(PETSC_COMM_WORLD, &x);
  VecSetSizes(x, PETSC_DECIDE, n);
  VecSetFromOptions(x);
  VecDuplicate(x, &y);

  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, n, n);
  MatSetFromOptions(A);
  MatSetUp(A);

  // Set A = I (identity)
  for (PetscInt i = 0; i < n; ++i) {
    PetscScalar one = 1.0;
    MatSetValue(A, i, i, one, INSERT_VALUES);
  }
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  VecSet(x, 2.0);
  MatMult(A, x, y);

  PetscReal norm = 0.0;
  VecNorm(y, NORM_2, &norm);

  int rank = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "  PETSc MatMult OK, ||y||_2 = " << norm << "\n";
  }

  MatDestroy(&A);
  VecDestroy(&x);
  VecDestroy(&y);
}
#endif

#ifdef FTDQMC_USE_SLEPC
static void test_slepc_eps_minimal() {
  EPS eps;
  EPSCreate(PETSC_COMM_WORLD, &eps);
  EPSSetProblemType(eps, EPS_NHEP);
  EPSDestroy(&eps);

  int rank = 0;
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  if (rank == 0) {
    std::cout << "  SLEPc EPSCreate/EPSDestroy OK\n";
  }
}
#endif

int main(int argc, char** argv) {
  // Initialization strategy:
  // - If SLEPc is enabled: SlepcInitialize() handles PETSc (+MPI) init.
  // - Else if PETSc is enabled: PetscInitialize() handles PETSc (+MPI) init.
  // - Else if MPI is enabled: MPI_Init().
  // This avoids double-initialization.

#ifdef FTDQMC_USE_MPI
  MPI_Init(&argc, &argv);
#endif

#ifdef FTDQMC_USE_SLEPC
  SlepcInitialize(&argc, &argv, nullptr, nullptr);
#elif defined(FTDQMC_USE_PETSC)
  PetscInitialize(&argc, &argv, nullptr, nullptr);
#endif

  int rank = 0;
#ifdef FTDQMC_USE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#elif defined(FTDQMC_USE_PETSC) || defined(FTDQMC_USE_SLEPC)
  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
#endif

  if (rank == 0) {
    print_build_flags();
  }

#ifdef FTDQMC_USE_OPENMP
  if (rank == 0) {
    std::cout << "  OpenMP max threads = " << omp_get_max_threads() << "\n";
  }
#endif

  // Always test BLAS
  test_blas_dgemm();

#ifdef FTDQMC_USE_MPI
  // If PETSc/SLEPc enabled, MPI is already initialized, but this still works
  test_mpi_allreduce();
#endif

#ifdef FTDQMC_USE_PETSC
  test_petsc_vec_mat();
#endif

#ifdef FTDQMC_USE_SLEPC
  test_slepc_eps_minimal();
#endif

  // Finalize
#ifdef FTDQMC_USE_SLEPC
  SlepcFinalize();
#elif defined(FTDQMC_USE_PETSC)
  PetscFinalize();
#elif defined(FTDQMC_USE_MPI)
  MPI_Finalize();
#endif

  return 0;
}