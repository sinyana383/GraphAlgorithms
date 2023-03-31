#include <gtest/gtest.h>

#include "../navigator/s21_graph/s21_graph.h"
#include "../navigator/s21_graph/s21_graph_algorithms.h"

class GraphAlgorithmsFixture : public ::testing::Test {
 protected:
  virtual void TearDown() override { delete _g0; }
  void SetUp() override {
    _g0 = new Graph();
    _g1 = new Graph();
  }

  Graph *_g0;
  Graph *_g1;
};

void check(std::vector<int> returned, std::vector<int> expected) {
  int i = 0;
  int size = returned.size();
  while (i < size) {
    EXPECT_EQ(returned[i], expected[i]);
    i += 1;
  }
}

TEST_F(GraphAlgorithmsFixture, depthFirstSearch) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/complex2.txt"));
  std::vector<int> arr = GraphAlgorithms::depthFirstSearch(*_g0, 1);
  char exp[] = {'a', 'd', 'c', 'b', 'e', 96, 96, 96, 96};
  char exp2[] = {'i', 'f', 96, 96, 96, 96, 96, 96, 96};
  char exp3[] = {'h', 96, 96, 96, 96, 96, 96, 96, 96};

  for (std::size_t i = 0; i < _g0->getVerticesNumber(); ++i)
    EXPECT_EQ(arr[i] + 'a' - 1, exp[i]);

  arr = GraphAlgorithms::depthFirstSearch(*_g0, 9);
  for (std::size_t i = 0; i < _g0->getVerticesNumber(); ++i)
    EXPECT_EQ(arr[i] + 'a' - 1, exp2[i]);

  arr = GraphAlgorithms::depthFirstSearch(*_g0, 8);
  for (std::size_t i = 0; i < _g0->getVerticesNumber(); ++i)
    EXPECT_EQ(arr[i] + 'a' - 1, exp3[i]);
}

TEST_F(GraphAlgorithmsFixture, getShortestPathBetweenVertices) {
  EXPECT_NO_THROW(_g1->loadGraphFromFile("unit_test_files/deicstra.txt"));
  EXPECT_EQ(60, GraphAlgorithms::getShortestPathBetweenVertices(*_g1, 1, 5));
  EXPECT_EQ(30, GraphAlgorithms::getShortestPathBetweenVertices(*_g1, 4, 5));

  EXPECT_EQ(std::numeric_limits<double>::infinity(),
            GraphAlgorithms::getShortestPathBetweenVertices(*_g1, 1, 6));

  EXPECT_EQ(0, GraphAlgorithms::getShortestPathBetweenVertices(*_g1, 1, 1));
}

TEST_F(GraphAlgorithmsFixture, getShortestPathsBetweenAllVertices) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/floid.txt"));
  EXPECT_NO_THROW(_g1->loadGraphFromFile("unit_test_files/floid_expected.txt"));

  std::vector<std::vector<double>> res =
      GraphAlgorithms::getShortestPathsBetweenAllVertices(*_g0);
  ASSERT_EQ(_g1->getVerticesNumber(), res.size());
  for (int i = 0; i < res.size(); ++i) {
    for (int j = 0; j < res[i].size(); ++j) {
      if (res[i][j] == std::numeric_limits<double>::infinity())
        EXPECT_EQ(0, _g1->getDist(i, j));
      else
        EXPECT_EQ(res[i][j], _g1->getDist(i, j));
    }
  }

  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/complex2.txt"));
  EXPECT_NO_THROW(
      _g1->loadGraphFromFile("unit_test_files/complex2_floid_exp.txt"));

  res = GraphAlgorithms::getShortestPathsBetweenAllVertices(*_g0);
  ASSERT_EQ(_g1->getVerticesNumber(), res.size());
  for (int i = 0; i < res.size(); ++i) {
    for (int j = 0; j < res[i].size(); ++j) {
      if (res[i][j] == std::numeric_limits<double>::infinity())
        EXPECT_EQ(0, _g1->getDist(i, j));
      else
        EXPECT_EQ(res[i][j], _g1->getDist(i, j));
    }
  }
}

TEST_F(GraphAlgorithmsFixture, breadthFirstSearch_1) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/complex.txt"));
  std::vector<int> arr = GraphAlgorithms::breadthFirstSearch(*_g0, 1);
  std::vector<int> returned;
  std::vector<int> expected = {1, 4, 3, 5, 2};

  returned = GraphAlgorithms::breadthFirstSearch(*_g0, 1);

  check(returned, expected);
}

TEST_F(GraphAlgorithmsFixture, breadthFirstSearch_2) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/complex.txt"));
  std::vector<int> arr = GraphAlgorithms::breadthFirstSearch(*_g0, 2);
  std::vector<int> returned;
  std::vector<int> expected = {2, 1, 5, 4, 3};

  returned = GraphAlgorithms::breadthFirstSearch(*_g0, 2);

  check(returned, expected);
}

TEST_F(GraphAlgorithmsFixture, breadthFirstSearch_3) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/complex.txt"));
  std::vector<int> arr = GraphAlgorithms::breadthFirstSearch(*_g0, 9);
  std::vector<int> returned;
  std::vector<int> expected = {9, 6, 7};

  returned = GraphAlgorithms::breadthFirstSearch(*_g0, 9);

  check(returned, expected);
}

TEST_F(GraphAlgorithmsFixture, TSM_test_nVerts_from_result_struct) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/saleman21.txt"));
  size_t nVert = _g0->getVerticesNumber() + 1;
  TsmResult res = GraphAlgorithms::solveTravelingSalesmanProblem(*_g0);
  EXPECT_EQ(nVert, res.vertices.size());
}

TEST_F(GraphAlgorithmsFixture, TSM_test_result_distance_gt_zero) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/saleman21.txt"));
  TsmResult res = GraphAlgorithms::solveTravelingSalesmanProblem(*_g0);
  EXPECT_GT(res.distance, 0);
}

TEST_F(GraphAlgorithmsFixture, TSM_test_errorFlag) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/salemanError.txt"));

  std::ofstream out("out.txt");
  std::streambuf *coutbuf = std::cout.rdbuf();  // save old buf
  std::cout.rdbuf(out.rdbuf());  // redirect std::cout to out.txt!

  TsmResult res = GraphAlgorithms::solveTravelingSalesmanProblem(*_g0);
  std::string test;
  std::ifstream f("out.txt");
  std::stringstream ss;
  ss << f.rdbuf();
  test = ss.str();
  std::cout.rdbuf(coutbuf);  // reset to standard output again
  remove("out.txt");
  EXPECT_EQ(test, "Error Graph\n");
}

TEST_F(GraphAlgorithmsFixture, MinOstovTree_check_resultMatrix) {
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/spanningTree7.txt"));

  std::vector<std::vector<int>> test;
  test.resize(_g0->getVerticesNumber());
  test[0] = {0, 0, 0, 0, 5, 1, 0};
  test[1] = {0, 0, 0, 0, 0, 0, 0};
  test[2] = {0, 3, 0, 0, 0, 0, 0};
  test[3] = {0, 0, 0, 0, 0, 0, 0};
  test[4] = {0, 0, 0, 0, 0, 0, 1};
  test[5] = {0, 0, 1, 1, 0, 0, 0};
  test[6] = {0, 0, 0, 0, 0, 0, 0};

  EXPECT_EQ(test, GraphAlgorithms::getLeastSpanningTree(*_g0));

  //-----------------------------------------------------------
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/spanningTree5.txt"));

  std::vector<std::vector<int>> test2;
  test2.resize(_g0->getVerticesNumber());
  test2[0] = {0, 0, 0, 2, 5};
  test2[1] = {0, 0, 0, 0, 0};
  test2[2] = {0, 3, 0, 0, 0};
  test2[3] = {0, 0, 1, 0, 0};
  test2[4] = {0, 0, 0, 0, 0};

  EXPECT_EQ(test2, GraphAlgorithms::getLeastSpanningTree(*_g0));

  //-----------------------------------------------------------
  EXPECT_NO_THROW(_g0->loadGraphFromFile("unit_test_files/spanningTree4.txt"));

  std::vector<std::vector<int>> test3;
  test3.resize(_g0->getVerticesNumber());
  test3[0] = {0, 0, 0, 2};
  test3[1] = {0, 0, 0, 0};
  test3[2] = {0, 0, 0, 0};
  test3[3] = {0, 3, 7, 0};

  EXPECT_EQ(test3, GraphAlgorithms::getLeastSpanningTree(*_g0));
}