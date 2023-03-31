#include <iostream>

#include "../navigator/s21_graph/s21_graph.h"
#include "../navigator/s21_graph/s21_graph_algorithms.h"

class Console_interface {
 private:
  Graph graph;

  auto Loading_file() -> void;
  auto read_file() -> void;
  auto to_do_list() -> void;
  auto depthFirstSearch() -> void;
  auto breadthFirstSearch() -> void;
  auto getShortestPathBetweenVertices() -> void;
  auto getShortestPathsBetweenAllVertices() -> void;
  auto getLeastSpanningTree() -> void;
  auto solveTravelingSalesmanProblem() -> void;
  auto checkVertice(int Vertex) -> bool;
  auto checkInput(int input, std::string str) -> int;

 public:
  Console_interface() {}
  ~Console_interface() {}

  auto Header_menu() -> void;
};
