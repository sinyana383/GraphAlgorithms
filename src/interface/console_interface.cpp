#include "console_interface.h"

#include <cstdlib>
#include <exception>
#include <iostream>

#include "../navigator/s21_graph/s21_graph_algorithms.h"

// using namespace std;

auto Console_interface::Header_menu() -> void {
  std::cout << "---------- WELCOME TO SIMPLE NAVIGATOR ----------" << std::endl;
  Loading_file();
}

auto Console_interface::Loading_file() -> void {
  std::cout << "ENTER FILE'S NAME: ";
  read_file();
}

auto Console_interface::read_file() -> void {
  std::string filename;
  std::cin >> filename;
  while (true) {
    try {
      graph.loadGraphFromFile(filename);
      break;
    } catch (std::exception& exc) {
      std::cerr << "\nError: " << exc.what() << "\n";
      Loading_file();
    }
  }
  to_do_list();
}

auto Console_interface::to_do_list() -> void {
  std::string str;
  int solve;

  while (1) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << " 1. loading the original graph from a file     " << std::endl;
    std::cout << " 2. breadthFirstSearch                         " << std::endl;
    std::cout << " 3. depthFirstSearch                           " << std::endl;
    std::cout << " 4. getShortestPathBetweenVertices             " << std::endl;
    std::cout << " 5. getShortestPathsBetweenAllVertices         " << std::endl;
    std::cout << " 6. getLeastSpanningTree                       " << std::endl;
    std::cout << " 7. solveTravelingSalesmanProblem              " << std::endl;
    std::cout << " 0. exit                                       " << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;

    solve = checkInput(solve, str);
    if (solve != -1) break;
  }
  switch (solve) {
    case 1:
      std::cout << "ENTER FILE'S NAME: ";
      read_file();
    case 2:
      breadthFirstSearch();
      to_do_list();
    case 3:
      depthFirstSearch();
      to_do_list();
    case 4:
      getShortestPathBetweenVertices();
      to_do_list();
    case 5:
      getShortestPathsBetweenAllVertices();
      to_do_list();
    case 6:
      getLeastSpanningTree();
      to_do_list();
    case 7:
      solveTravelingSalesmanProblem();
      to_do_list();
    case 0:
      std::cout << "See you! Bye bye.\n" << std::endl;
      exit(0);
      break;
    default: {
      std::cerr << "Error: incorrect input, try again " << std::endl;
      to_do_list();
    }
  }
}

auto Console_interface::depthFirstSearch() -> void {
  int vertex;
  std::string str;

  while (1) {
    std::cout << "Please, enter number of vertex\n" << std::endl;
    vertex = checkInput(vertex, str);
    if (vertex != -1 && checkVertice(vertex))
      break;
    else
      std::cout << "That's wrong, enter another number of vertex\n"
                << std::endl;
  }
  std::vector<int> arr = GraphAlgorithms::depthFirstSearch(graph, vertex);
  for (auto i : arr) {
    if (i == 0) break;
    std::cout << i << " ";
  }
  std::cout << std::endl;
}

auto Console_interface::breadthFirstSearch() -> void {
  int vertex;
  std::string str;

  while (1) {
    std::cout << "Please, enter number of vertex\n" << std::endl;
    vertex = checkInput(vertex, str);
    if (vertex != -1 && checkVertice(vertex))
      break;
    else
      std::cout << "That's wrong, enter another number of vertex\n"
                << std::endl;
  }
  std::vector<int> arr = GraphAlgorithms::breadthFirstSearch(graph, vertex);
  for (auto i : arr) {
    std::cout << i << " ";
  }
  std::cout << std::endl;
}

auto Console_interface::getShortestPathBetweenVertices() -> void {
  int vertex_1;
  int vertex_2;
  std::string str;

  while (1) {
    std::cout << "Please, enter first vertex\n" << std::endl;
    vertex_1 = checkInput(vertex_1, str);
    if (vertex_1 != -1 && checkVertice(vertex_1))
      break;
    else
      std::cout << "That's wrong, enter another number of vertex\n"
                << std::endl;
  }
  while (1) {
    std::cout << "Please, enter second vertex\n" << std::endl;
    vertex_2 = checkInput(vertex_2, str);
    if (vertex_2 != -1 && checkVertice(vertex_2))
      break;
    else
      std::cout << "That's wrong, enter another number of vertex\n"
                << std::endl;
  }
  std::cout << "the shortest path is "
            << GraphAlgorithms::getShortestPathBetweenVertices(graph, vertex_1,
                                                               vertex_2)
            << std::endl;
}

auto Console_interface::getShortestPathsBetweenAllVertices() -> void {
  std::cout << "The Shortest Paths Between All Vertices:\n" << std::endl;
  std::vector<std::vector<double>> arr =
      GraphAlgorithms::getShortestPathsBetweenAllVertices(graph);
  for (auto i : arr) {
    for (auto el : i) {
      if (el == std::numeric_limits<double>::infinity())
        std::cout << 0 << " ";
      else
        std::cout << el << " ";
    }
    std::cout << std::endl;
  }
}

auto Console_interface::getLeastSpanningTree() -> void {
  std::cout << "Spanning Tree\n" << std::endl;
  std::vector<std::vector<int>> mtrx;
  mtrx = GraphAlgorithms::getLeastSpanningTree(graph);

  {
    size_t size = graph.getVerticesNumber();

    for (size_t i = 0; i < size; i++) {
      for (size_t j = 0; j < size; j++) {
        printf("%d ", mtrx[i][j]);
      }
      printf("\n");
    }
  }
}

auto Console_interface::solveTravelingSalesmanProblem() -> void {
  extern bool g_errorGraphTsm;
  TsmResult tsm;

  tsm = GraphAlgorithms::solveTravelingSalesmanProblem(graph);
  if (g_errorGraphTsm) return;
  std::cout << "Solve Salesman's problem is\n" << std::endl;

  {
    printf("Way:");
    for (size_t i = 0; i < tsm.vertices.size(); i++) {
      printf(" %d", tsm.vertices[i]);
    }
    printf("\n");
    printf("Distance= %f\n", tsm.distance);
  }
}

auto Console_interface::checkVertice(int Vertex) -> bool {
  if (Vertex < 1 || Vertex > graph.getVerticesNumber()) return 0;
  return 1;
}

auto Console_interface::checkInput(int input, std::string str) -> int {
  std::size_t pos{};
  std::cin >> str;
  try {
    input = std::stoi(str, &pos);
  } catch (std::invalid_argument const& ex) {
    std::cout << "std::invalid_argument::what(): " << ex.what() << '\n';
    std::cin.clear();
    return (-1);
  } catch (std::out_of_range const& ex) {
    std::cout << "std::out_of_range::what(): " << ex.what() << '\n';
    std::cin.clear();
    return (-1);
  }
  if (std::cin.peek() != '\n') {
    std::cerr << "Error: incorrect input, try again " << std::endl;
    std::cin.clear();
    return (-1);
  }
  std::cin.clear();
  return input;
}