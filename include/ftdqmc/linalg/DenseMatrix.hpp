#pragma once

#include <cmath>
#include <complex>
#include <cstddef>
#include <vector>

namespace ftdqmc::linalg {

using Complex = std::complex<double>;

class DenseMatrix {
 public:
  DenseMatrix() = default;
  DenseMatrix(int rows, int cols)
      : rows_(rows), cols_(cols), ld_(rows), data_(static_cast<size_t>(ld_) * cols) {}

  int rows() const { return rows_; }
  int cols() const { return cols_; }
  int ld() const { return ld_; }

  Complex& operator()(int i, int j) { return data_[index(i, j)]; }
  const Complex& operator()(int i, int j) const { return data_[index(i, j)]; }

  static DenseMatrix zeros(int rows, int cols) { return DenseMatrix(rows, cols); }

 private:
  int rows_ = 0;
  int cols_ = 0;
  int ld_ = 0;
  std::vector<Complex> data_; // column-major: index = i + ld * j

  size_t index(int i, int j) const {
    return static_cast<size_t>(i) + static_cast<size_t>(ld_) * static_cast<size_t>(j);
  }
};

inline Complex trace(const DenseMatrix& m) {
  const int n = (m.rows() < m.cols()) ? m.rows() : m.cols();
  Complex tr{0.0, 0.0};
  for (int i = 0; i < n; ++i) {
    tr += m(i, i);
  }
  return tr;
}

inline double hermitian_error_norm(const DenseMatrix& m) {
  double sum = 0.0;
  for (int i = 0; i < m.rows(); ++i) {
    for (int j = 0; j < m.cols(); ++j) {
      const Complex diff = m(i, j) - std::conj(m(j, i));
      sum += std::norm(diff);
    }
  }
  return std::sqrt(sum);
}

}  // namespace ftdqmc::linalg
