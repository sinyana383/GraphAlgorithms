#include "Matrix.hpp"
Matrix::Matrix(std::size_t size) : _sideSize(size) {
  _matrix = new double*[_sideSize];

  for (int i = 0; i < _sideSize; ++i) _matrix[i] = new double[_sideSize];
}
Matrix::~Matrix() {
  for (int i = 0; i < _sideSize; ++i) delete[] _matrix[i];
  delete[] _matrix;
}
