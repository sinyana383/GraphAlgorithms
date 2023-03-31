#ifndef NAVIGATOR_MATRIX_HPP
#define NAVIGATOR_MATRIX_HPP

#include <cstdio>
class Matrix {
 private:
  std::size_t _sideSize;
  double** _matrix;

 public:
  explicit Matrix(std::size_t size);
  ~Matrix();

  std::size_t getSideSize() const { return _sideSize; }

  double& at(size_t row, size_t col) { return _matrix[row][col]; }
};

#endif  // NAVIGATOR_MATRIX_HPP
