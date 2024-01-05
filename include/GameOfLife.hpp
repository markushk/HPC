#ifndef GAMEOFLIFE_H
#define GAMEOFLIFE_H

#include <vector>

class GameOfLife {
public:
    GameOfLife(int rows, int cols, int seed, double probability);
    ~GameOfLife();

    void initialConfiguration();
    void printStatus();
    void printFieldAnimated();
    void printField();
    void runLife(int generations);

private:
    int _rows;
    int _cols;
    int _seed;
    double _probability;
    std::vector<std::vector<char>> _world;
    std::vector<std::vector<char>> _worldCopy;

    char getRandomValue(int row_id, int col_id);
    int countAlive();
    int countNeighbours(int i, int j);

    int plus(int x, int m);
    int minus(int x, int m);
};

#endif // GAMEOFLIFE_H
