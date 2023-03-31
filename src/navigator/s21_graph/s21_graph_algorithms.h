#ifndef NAVIGATOR_SRC_GRAPHALGORITHMS_H_
#define NAVIGATOR_SRC_GRAPHALGORITHMS_H_

#include <limits>
#include <vector>

#include "../Queue/s21_queue.hpp"
#include "../stack/Stack.hpp"
#include "s21_graph.h"

struct TsmResult {
  std::vector<int> vertices;  // an array with the route you are looking for
                              // (with the vertex traverse order). Instead of
                              // int* you can use std::vector<int>
  double distance;            // the length of this route
};

class GraphAlgorithms {
 public:
  static std::vector<int> breadthFirstSearch(Graph &graph,
                                             int startVertex);              // 2
  static std::vector<int> depthFirstSearch(Graph &graph, int startVertex);  // 1

  static double getShortestPathBetweenVertices(Graph &graph, int vertex1,
                                               int vertex2);  // 3
  static std::vector<std::vector<double> > getShortestPathsBetweenAllVertices(
      Graph &graph);  // 4
  static std::vector<std::vector<int> > getLeastSpanningTree(
      Graph &graph);                                             // 5
  static TsmResult solveTravelingSalesmanProblem(Graph &graph);  // 6
};

#endif  // NAVIGATOR_SRC_GRAPHALGORITHMS_H_
