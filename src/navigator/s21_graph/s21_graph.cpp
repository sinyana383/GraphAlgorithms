#include "s21_graph.h"

void Graph::loadGraphFromFile(const std::string &filename) {
  std::ifstream fin;
  fin.open(filename, std::ios::in);
  if (!fin.is_open() || fin.bad()) throw std::runtime_error("loadFile error");

  std::string line;
  std::getline(fin, line);
  std::size_t size = std::stoi(line);
  delete _m;
  _type = UNDIRECT;
  _m = new Matrix(size);

  char *pStart;
  std::size_t pEnd;
  int inCol, inRow = -1;
  while (std::getline(fin, line) && ++inRow < size) {
    inCol = 0;
    pStart = &(line[0]);
    do {
      try {
        _m->at(inRow, inCol) = std::stoi(pStart, &pEnd);
      } catch (std::exception &e) {
        throw std::runtime_error(e.what());
      }
      if (_m->at(inRow, inCol) < 0)
        throw std::runtime_error("negative distance");
      if (inCol < inRow && _m->at(inCol, inRow) != _m->at(inRow, inCol))
        _type = DIRECT;
      pStart += pEnd;
    } while (*pStart && ++inCol < size);
    if (inCol != size - 1) throw std::runtime_error("wrong columns number");
  }
  if (inRow != size - 1) throw std::runtime_error("wrong rows number");
  fin.close();
}
int Graph::exportGraphToDot(const std::string &filename) {
  std::ofstream fout;
  fout.open(filename, std::ios::out | std::ios::trunc);
  if (!fout.is_open() || fout.bad()) return error("exportFile error");

  if (_type == UNDIRECT)
    fout << "graph";
  else
    fout << "digraph";
  fout << " graphname {" << std::endl;
  for (int i = 0; i < _m->getSideSize(); ++i)
    fout << "    " << (char)('a' + i) << ";" << std::endl;

  for (int r = 0; r < _m->getSideSize(); ++r) {
    for (int c = 0; c < _m->getSideSize(); ++c) {
      if (_m->at(r, c) > 0) {
        if (_m->at(c, r) == _m->at(r, c) && r < c)
          fout << "    " << (char)('a' + r) << " -- " << (char)('a' + c) << ";"
               << std::endl;
        else if (_m->at(c, r) != _m->at(r, c))
          fout << "    " << (char)('a' + r) << " -> " << (char)('a' + c) << ";"
               << std::endl;
      }
    }
  }
  fout << "}";
  fout.close();
  return 1;
}

double Graph::getDist(size_t a, size_t b) {
  if (a >= _m->getSideSize() || b >= _m->getSideSize())
    throw std::out_of_range("wrong vertices number");
  return _m->at(a, b);
}

int error(const std::string &massage) {
  std::cout << massage << std::endl;
  return -1;
}
