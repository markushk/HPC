#include "GameOfLife.hpp"
#include <iostream>
#include <cstdlib>
#include <thread>
#include <chrono>

GameOfLife::GameOfLife(int rows, int cols, int seed, double probability)
    : _rows(rows), _cols(cols), _seed(seed), _probability(probability),
      _world(rows, std::vector<char>(cols)), _worldCopy(rows, std::vector<char>(cols)) {}

GameOfLife::~GameOfLife() {}

void GameOfLife::initialConfiguration() {
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            _world[i][j] = getRandomValue(i, j);
        }
    }
}

void GameOfLife::printStatus() {
    int alive = countAlive();
    std::cout << "Number of alive cells: " << alive << std::endl;
    std::cout << "Number of dead cells: " << _rows * _cols - alive << std::endl;
}

void GameOfLife::printFieldAnimated() {
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            if (_world[i][j] == '1')
                std::cout << " X";
            else
                std::cout << " .";
        }
        std::cout << "\n";
    }
    std::cout << "\x1B[" << _rows << "A";  // Carriage return to move back to the beginning of the line
    std::this_thread::sleep_for(std::chrono::milliseconds(300));  // Adjust delay as needed
    std::cout.flush();  // Flush the output
}

void GameOfLife::printField() {
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            if (_world[i][j] == '1')
                std::cout << " X";
            else
                std::cout << " .";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void GameOfLife::runLife(int generations) {
    for (int gen = 0; gen < generations; gen++) {
        for (int i = 0; i < _rows; ++i) {
            for (int j = 0; j < _cols; ++j) {
                int neighbours = countNeighbours(i, j);
                if (_world[i][j] == '0') {
                    if (neighbours == 3) {
                        _worldCopy[i][j] = '1';
                    } else {
                        _worldCopy[i][j] = '0';
                    }
                } else {
                    if (neighbours != 2 && neighbours != 3) {
                        _worldCopy[i][j] = '0';
                    } else {
                        _worldCopy[i][j] = '1';
                    }
                }
            }
        }
        _world = _worldCopy;
    }
}

char GameOfLife::getRandomValue(int row_id, int col_id) {
    char r = '0';
    int my_seed = _seed + row_id * _cols + col_id;
    srand(my_seed);

    double rand_value = static_cast<double>(rand()) / RAND_MAX;

    if (rand_value < _probability) {
        r = '1';
    }

    return r;
}

int GameOfLife::countAlive() {
    int alive = 0;
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            if (_world[i][j] == '1')
                alive++;
        }
    }
    return alive;
}

int GameOfLife::plus(int x, int m) {
    return (x + 1) % m;
}

int GameOfLife::minus(int x, int m) {
    return (x - 1 + m) % m;
}

int GameOfLife::countNeighbours(int i, int j) {
    int number_of_neighbours = 0;
    if (_world[i][plus(j, _cols)] == '1') number_of_neighbours++;     // right
    if (_world[minus(i, _rows)][plus(j, _cols)] == '1') number_of_neighbours++;  // down right
    if (_world[minus(i, _rows)][j] == '1') number_of_neighbours++;    // down
    if (_world[minus(i, _rows)][minus(j, _cols)] == '1') number_of_neighbours++;  // down left
    if (_world[i][minus(j, _cols)] == '1') number_of_neighbours++;    // left
    if (_world[plus(i, _rows)][minus(j, _cols)] == '1') number_of_neighbours++;   // up left
    if (_world[plus(i, _rows)][j] == '1') number_of_neighbours++;     // up
    if (_world[plus(i, _rows)][plus(j, _cols)] == '1') number_of_neighbours++;  // up right
    return number_of_neighbours;
}


