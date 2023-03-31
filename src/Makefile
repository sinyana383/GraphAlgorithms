CC = g++

CFLAGS = -Wall -Werror -Wextra --std=c++17
BUILD_DIR = build
CLANG_TIDY_CMD = clang-format -style=google -n

.PHONY: all clean test install lint

all: s21_graph.a s21_graph_algorithms.a install

$(BUILD_DIR):
	mkdir -p $@

s21_graph.a:
	cmake -S . -B build && cmake --build build --target s21_graph
	mv build/libs21_graph.a ./s21_graph.a

s21_graph_algorithms.a:
	cmake -S . -B build && cmake --build build --target s21_graph_algorithms
	mv build/libs21_graph_algorithms.a ./s21_graph_algorithms.a

install: $(BUILD_DIR)
	cmake -S . -B build && cmake --build build --target navigator
	./$(BUILD_DIR)/navigator

lint:
	$(CLANG_TIDY_CMD) navigator/*/* interface/* test/*.cpp navigator/main.cpp

test:
	cmake -S . -B build && cmake --build build --target testing
	./$(BUILD_DIR)/testing

clean:
	rm -rf $(BUILD_DIR) s21_graph.a s21_graph_algorithms.a