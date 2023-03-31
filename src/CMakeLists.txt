cmake_minimum_required(VERSION 3.18)
project(navigator)

set(CMAKE_CXX_STANDARD 17)

####LIBRARIES####

# make s21_graph.a
add_library(s21_graph STATIC
        # файлы по s21_graph
        navigator/s21_graph/s21_graph.h 
        navigator/s21_graph/s21_graph.cpp
        navigator/matrix/Matrix.cpp navigator/matrix/Matrix.hpp
        )

# make s21_graph_algorithms.a
add_library(s21_graph_algorithms STATIC 
        # файлы по s21_graph_algorithms
        navigator/s21_graph/s21_graph_algorithms.h 
        navigator/s21_graph/s21_graph_algorithms.cpp
        navigator/stack/Stack.hpp navigator/Queue/s21_queue.hpp)

####EXECUTABLES####

# make test
add_subdirectory(test/googletest-main)
include_directories(test/googletest-main/googletest/include)
add_executable(testing
        # файлы по unit-test
        test/graphTest.cpp
        test/graphAlgorithmsTest.cpp
        )
target_link_libraries(testing s21_graph_algorithms s21_graph gtest)

# make install
include_directories(navigator/s21_graph)
add_executable(navigator
        # файлы по интерфейсу
        navigator/main.cpp
        interface/console_interface.h
        interface/console_interface.cpp
)
target_link_libraries(navigator s21_graph_algorithms s21_graph)
