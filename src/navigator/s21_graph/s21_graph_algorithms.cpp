#include "s21_graph_algorithms.h"

#include <algorithm>
#include <cassert>
#include <cmath>
bool g_errorGraphTsm;

std::vector<int> GraphAlgorithms::depthFirstSearch(Graph &graph,
                                                   int startVertex) {
  int iRes = 0;
  std::vector<int> res(graph.getVerticesNumber(), 0);
  Stack<int> stack;
  std::vector<int> visited(graph.getVerticesNumber(), 0);
  stack.init();

  stack.push(startVertex - 1);
  visited[stack.peek()] = 1;  // защита от зацикленных пунктов
  res[iRes++] = stack.peek() + 1;
  while (stack.GetSize() > 0) {
    for (int c = 0; c < graph.getVerticesNumber(); ++c) {
      if (graph.getDist(stack.peek(), c) != 0 && visited[c] == 0) {
        visited[c] = 1;
        res[iRes++] = c + 1;
        stack.push(c);
        break;
      }
      if (c == graph.getVerticesNumber() - 1) stack.pop();
    }
  }
  return res;
}

double GraphAlgorithms::getShortestPathBetweenVertices(Graph &graph,
                                                       int vertex1,
                                                       int vertex2) {
  std::vector<double> paths(graph.getVerticesNumber(),
                            std::numeric_limits<double>::infinity());
  std::vector<bool> visited(graph.getVerticesNumber(), false);

  auto minEl = paths.begin() + vertex1 - 1;
  int cur = vertex1 - 1;
  paths[cur] = 0;
  visited[cur] = true;

  // защита и на vertex1 == vertex2 и на отсутствие пути к вершине
  for (int finds = 0; finds < graph.getVerticesNumber() && cur != vertex2 - 1;
       ++finds) {
    for (int c = 0; c < graph.getVerticesNumber(); ++c) {
      double dist = graph.getDist(cur, c);
      if (!(dist == 0 || cur == c || visited[c]))
        paths[c] = std::min(paths[c], paths[cur] + dist);

      if ((paths[c] - *minEl < 0 || *minEl == 0) && !visited[c])
        minEl = paths.begin() + c;
    }
    cur = minEl - paths.begin();
    minEl = paths.begin() + vertex1 - 1;
    visited[cur] = true;
  }
  if (cur == vertex2 - 1) return paths[cur];
  return std::numeric_limits<double>::infinity();
}
std::vector<std::vector<double>>
GraphAlgorithms::getShortestPathsBetweenAllVertices(Graph &graph) {
  std::size_t sideSize = graph.getVerticesNumber();
  std::vector<std::vector<double>> memo(sideSize);
  std::vector<std::vector<double>> prev(sideSize);
  for (int i = 0; i < sideSize; ++i) {
    memo[i].resize(sideSize);
    prev[i].resize(sideSize);
    for (int j = 0; j < sideSize; ++j) {
      if (i == j)
        prev[i][j] = 0;
      else if (graph.getDist(i, j) == 0)
        prev[i][j] = std::numeric_limits<double>::infinity();
      else
        prev[i][j] = graph.getDist(i, j);
    }
  }

  for (int k = 0; k < sideSize; ++k) {
    for (int i = 0; i < sideSize; ++i) {
      for (int j = 0; j < sideSize; ++j) {
        memo[i][j] = std::min(prev[i][j], prev[i][k] + prev[k][j]);
        prev[i][j] = memo[i][j];
      }
    }
  }
  return memo;
}

std::vector<int> GraphAlgorithms::breadthFirstSearch(
    Graph &graph,
    int startVertex) {  // areel

  std::vector<int> visited;
  int sizeG = graph.getVerticesNumber();
  ;  // sizeG - количество вершин графа
  int vertex, b;
  //	массив размера sizeG bool* passed // посещенные вершины
  bool *passed = new bool[sizeG];
  Queue<int> queue;
  queue.init();  // создание пустого списка очереди
  queue.push(startVertex - 1);  // добавляем заданную вершину в очередь, -1 т.к.
                                // начинаются с 1 в графе, а в с++ с 0

  while (queue.isEmpty() == false)  // пока в очереди есть вершины
  {
    vertex = queue.pop();  // получение эл-та из очереди с удалением -- это
                           // строка а в матрице смежности
    // если вершина уже посещена, то ничего не делаем, continue;
    if (passed[vertex] == true) {
      continue;
    }
    // делаем для вершины passed true;
    passed[vertex] = true;
    visited.push_back(vertex + 1);  // вставить новый элемент в конец вектора,
                                    // который увеличивает размер вектора на 1
    // graph.getDist(a,b) -  массив смежности используем graph.getDist
    b = 0;
    while (b < sizeG) {
      // если не посещена i вершина && смежные graph.getDist(queue.peek(),b) !=
      // 0
      if (passed[b] != true && graph.getDist(vertex, b) != 0) {
        queue.push(b);  // добавляем вершину в очередь
      }
      b++;
    }
  }
  delete[] passed;
  return visited;
}

//------------MinOstovTree:-----------------// Женя
template <typename T>
void remove(std::vector<T> &v, size_t index) {
  v.erase(v.begin() + index);
}

int findStartVert(std::vector<std::vector<int>> &helpMatrix, size_t nVertices) {
  for (size_t row = 0; row < nVertices; row++) {
    for (size_t col = 0; col < nVertices; col++) {
      if (helpMatrix[row][col] != 0) return row;
    }
  }
  return (-1);
}

bool checkAllElemsAreZero(std::vector<int> v) {
  for (size_t i = 0; i < v.size(); i++) {
    if (v[i] != 0) return (0);
  }
  return (1);
}

// Алгоритм Прима — алгоритм построения минимального остовного дерева
// взвешенного связного неориентированного графа.
std::vector<std::vector<int>> GraphAlgorithms::getLeastSpanningTree(
    Graph &graph) {
  size_t nVertices = graph.getVerticesNumber();
  std::vector<std::vector<int>> mtrxOstovTree(nVertices);
  std::vector<std::vector<int>> helpMatrix(nVertices);

  for (size_t i = 0; i < nVertices; i++) {
    mtrxOstovTree[i].resize(nVertices);
    helpMatrix[i].resize(nVertices);
  }

  for (size_t i = 0; i < nVertices; i++) {
    for (size_t j = 0; j < nVertices; j++) {
      // Избавляюсь от петлей (зануляю диагональные элементы)
      // заполняю вспомогательную матрицу (копия загруженной из файла)
      // заполняю нулями матрицу остовного дерева
      if (i != j)
        helpMatrix[i][j] = graph.getDist(i, j);
      else {
        helpMatrix[i][j] = 0;
      }
      mtrxOstovTree[i][j] = 0;
    }
  }
  /*	vector для хранения  вершин. На каждом шаге буду добавлять
           вершину в этот вектор, при построении остова дерева
  */
  std::vector<int> connectedVert(0);
  // ищем первую ненулевую вершину
  int startVert = findStartVert(helpMatrix, nVertices);
  connectedVert.push_back(startVert);

  // Search ostovTree from the 1st vertex
  // Stop search when all vertices will be connected;
  // на каждой итерации цикла в матрицу остовного дерева добавляется одна
  // вершина
  int tempRow, tempCol, iTemp;
  int n = 0;
  int minWeight;
  while (n < nVertices - 1) {
    minWeight = std::numeric_limits<int>::max();
    for (size_t i = 0; i < connectedVert.size(); i++) {
      int row = connectedVert.at(i);
      for (size_t col = 0; col < nVertices; col++) {
        if (std::find(connectedVert.begin(), connectedVert.end(), col) !=
            connectedVert.end())
          continue;
        if (helpMatrix[row][col] > 0 && helpMatrix[row][col] < minWeight) {
          minWeight = helpMatrix[row][col];
          tempRow = row;
          tempCol = col;
          iTemp = i;
        }
      }
    }
    mtrxOstovTree[tempRow][tempCol] = helpMatrix[tempRow][tempCol];
    helpMatrix[tempRow][tempCol] = 0;
    helpMatrix[tempCol][tempRow] = 0;

    if (checkAllElemsAreZero(helpMatrix[tempRow])) remove(connectedVert, iTemp);
    if (!checkAllElemsAreZero(helpMatrix[tempCol]))
      connectedVert.push_back(tempCol);  // можно оставить только эту строчку. а
                                         // эти две проверки удалить
    n++;
  }
  return (mtrxOstovTree);
}

//------------------------задача коммивояжера-------------

/*
        Был открыт ряд комбинаций α/β, которые позволяют находить хорошие
   результаты: α/β: 0.5/5   1/1   1/2   1/5; Rho (ρ): от 0 до 1
*/
typedef struct Data {
  const double alpha = 1;  // влияет на фермент
  const double beta = 1;   // на расстояние
  const double rho = 0.5;  // rho это распыление на грань; (1-rho) это испарение
  const double qVal = 1000;
  double initialPheromone;
  size_t nAnts;
  size_t nVerts;
  std::vector<std::vector<double>> pheromone;
  std::vector<int> bestWay;
  int bestDistance;
} Data;

typedef struct AntStruct {
  int curVert;
  int nextVert;
  std::vector<int> visited;
  std::vector<int> unvisited;
  double length;
} AntStruct;

bool errorGraph(Graph &graph, int nVerts) {
  for (size_t i = 0; i < nVerts; i++) {
    for (size_t n = 0; n < nVerts; n++) {
      if ((graph.getDist(i, n) == 0 && (i != n)) ||
          (graph.getDist(i, n) != graph.getDist(n, i))) {
        return 1;
      }
    }
  }
  return 0;
}

int getNextVert(AntStruct &ant, Graph &graph, Data &dataStruct) {
  ant.curVert = *(ant.visited.end() - 1);
  double denominator = 0.0;
  for (size_t unvisitInd = 0; unvisitInd < ant.unvisited.size();
       unvisitInd++) {  // для каждой непосещённой вершины считаю знаменатель:
    double distance = graph.getDist(ant.curVert, ant.unvisited[unvisitInd]);
    denominator +=
        pow(dataStruct.pheromone[ant.curVert][ant.unvisited[unvisitInd]],
            dataStruct.alpha) *
        pow(1.0 / distance, dataStruct.beta);
  }

  int unvisitInd = 0;
  while (1) {  // считаю числитель и вероятность; кидаю монетку , получаю
               // следующую вершину.
    double p, n;
    if (unvisitInd == ant.unvisited.size()) unvisitInd = 0;
    double distance = graph.getDist(ant.curVert, ant.unvisited[unvisitInd]);
    double nominator =
        pow(dataStruct.pheromone[ant.curVert][ant.unvisited[unvisitInd]],
            dataStruct.alpha) *
        pow(1.0 / distance, dataStruct.beta);
    if (denominator <= 0.0) denominator = 100 * nominator;
    assert(denominator != 0);
    p = 100 * nominator / denominator;
    if ((n = rand() % 100) < p) break;
    unvisitInd += 1;
  }
  int ret = ant.unvisited[unvisitInd];
  remove(ant.unvisited, unvisitInd);  // ant start from curVert to nextVert
  return ret;
}

int antsGoGoGo(AntStruct *ants, Graph &graph, Data &dataStruct) {
  // для каждого муравья:
  for (size_t i = 0; i < dataStruct.nAnts; i++) {
    for (; ants[i].unvisited.size() >
           0;) {  // на каждой итерации вектор unvisited уменьшается на одну
                  // вершину
      // выбрать след.вершину из вектора unvisited:
      ants[i].nextVert = getNextVert(ants[i], graph, dataStruct);
      //-----след.вершина выбрана
      ants[i].length += graph.getDist(ants[i].curVert, ants[i].nextVert);
      ants[i].visited.push_back(
          ants[i].nextVert);  // движение из curVert в nextVert
      // удаление из вектора ants[i].unvisited было в функции getNextVert()---
    }
    // добавить начальную вершину(ants[i].visited[0]):
    ants[i].curVert = *(ants[i].visited.end() - 1);
    ants[i].length += graph.getDist(ants[i].curVert, ants[i].visited[0]);
    ants[i].visited.push_back(ants[i].visited[0]);
    //---------------
    if (ants[i].length < dataStruct.bestDistance) {
      dataStruct.bestDistance = ants[i].length;
      dataStruct.bestWay.assign(ants[i].visited.begin(), ants[i].visited.end());
    }
  }
  return 1;
}

void restartAnts(AntStruct *ants, Data &dataStuct) {
  for (size_t i = 0; i < dataStuct.nVerts; i++) {
    ants[i].length = 0;
    ants[i].visited.clear();
    ants[i].unvisited.clear();
    ants[i].visited.push_back(i);  // муравей в каждой вершине
    for (size_t n = 0; n < dataStuct.nVerts; n++) {
      if (n == i) continue;
      ants[i].unvisited.push_back(n);
    }
  }
}

void getNewValues(AntStruct *ants, Data &dataStruct) {
  double deltaPherom, newPherom, oldPherom;

  // испарим фермент на каждой грани:
  for (size_t row = 0; row < dataStruct.nVerts; row++) {
    for (size_t col = 0; col < dataStruct.nVerts; col++) {
      if (col == row) continue;
      dataStruct.pheromone[row][col] =
          dataStruct.pheromone[row][col] * (1 - dataStruct.rho);
      if (dataStruct.pheromone[row][col] <= 0.0)
        dataStruct.pheromone[row][col] = dataStruct.initialPheromone;
    }
  }

  for (size_t i = 0; i < dataStruct.nAnts;
       i++) {  // для каждого муравья расчитать и разложить феромон по граням
               // пути муравья
    deltaPherom = (dataStruct.qVal /
                   ants[i].length);  // распылил i-муравей  dataStruct.rho *
    for (size_t n = 0; n < dataStruct.nVerts;
         n++) {  // для каждой грани пути текущего муравья:
      oldPherom = dataStruct.pheromone[*(ants[i].visited.begin() + n)]
                                      [*(ants[i].visited.begin() + n + 1)];
      newPherom = oldPherom + deltaPherom;
      newPherom = (newPherom <= 0.0) ? dataStruct.initialPheromone : newPherom;
      // записываем новый феромон в вектор-феромон:
      dataStruct.pheromone[*(ants[i].visited.begin() + n)]
                          [*(ants[i].visited.begin() + n + 1)] = newPherom;
    }
  }
}

TsmResult GraphAlgorithms::solveTravelingSalesmanProblem(Graph &graph) {
  TsmResult ret;
  size_t nVerts = graph.getVerticesNumber();
  g_errorGraphTsm = 0;
  if (errorGraph(graph, nVerts)) {
    std::cout << "Error Graph" << std::endl;
    g_errorGraphTsm = 1;
    return ret;
  }
  size_t nTimes = 20 * nVerts;
  AntStruct ants[nVerts];

  // Init:
  Data dataStruct;
  dataStruct.bestDistance = std::numeric_limits<int>::max();
  dataStruct.nAnts = nVerts;
  dataStruct.nVerts = nVerts;
  dataStruct.initialPheromone = 1.0 / (double)(nVerts);
  dataStruct.pheromone.resize(nVerts);
  for (size_t i = 0; i < nVerts; i++) {
    dataStruct.pheromone[i].resize(nVerts);
    for (size_t n = 0; n < nVerts; n++) {
      dataStruct.pheromone[i][n] = dataStruct.initialPheromone;
    }
  }

  for (size_t i = 0; i < nVerts; i++) {
    ants[i].visited.push_back(i);  // муравей в каждой вершине
    for (size_t n = 0; n < nVerts; n++) {
      if (n == i) continue;
      ants[i].unvisited.push_back(n);
    }
    ants[i].length = 0;
  }
  //--endInit---

  srand(time(NULL));
  for (size_t i = 0; i < nTimes; i++) {
    antsGoGoGo(ants, graph, dataStruct);
    getNewValues(ants, dataStruct);
    restartAnts(ants, dataStruct);
  }

  for (size_t i = 0; i < dataStruct.bestWay.size();
       i++) {  // по сабжекту, вершины должны начинаться не с 0, а с 1.
    dataStruct.bestWay[i] += 1;
  }
  ret.vertices.assign(dataStruct.bestWay.begin(), dataStruct.bestWay.end());
  ret.distance = dataStruct.bestDistance;
  return ret;
}
