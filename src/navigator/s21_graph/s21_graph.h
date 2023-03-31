#ifndef NAVIGATOR_S21_GRAPH_H
#define NAVIGATOR_S21_GRAPH_H

#include <fstream>
#include <iostream>
#include <string>

#include "../matrix/Matrix.hpp"

#define UNDIRECT 0
#define DIRECT 1

#define FILENAME "complex_test.txt"
#define DOTFILENAME "complex_test.dot"

class Graph {
 private:
  Matrix *_m;
  int _type = UNDIRECT;

 public:
  ~Graph() { delete _m; }
  Graph() : _m(new Matrix(0)) {}

  void loadGraphFromFile(
      const std::string &filename);  // загрузка матрицы смежности из файла
  int exportGraphToDot(const std::string &filename);  // выгрузка матрицы в фейл

  double getDist(size_t a,
                 size_t b);  // возвращает элемент матрицы: а/i - строка, b/j -
                             // столбец, нумерация с нуля
  std::size_t getVerticesNumber() {
    return _m->getSideSize();
  }  // кол-во вершин в графе(т. е. длина и ширина матрицы)
};

int error(const std::string &massage);

#endif  // NAVIGATOR_S21_GRAPH_H
