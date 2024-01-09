#ifndef GAMEOFLIFE_H
#define GAMEOFLIFE_H

#include <vector>
#include <mpi.h>

class GameOfLife {
public:
    GameOfLife(int rows, int cols, int seed, double probability);
    ~GameOfLife();

    void initialConfiguration();
    void printStatus();
    void printFieldAnimated();
    void printField();
    void runLife(int generations);
    void init_mpi();
    void finalizempi();
    void testCommunication();
    void pointToPoint();

private:
    int _rows;
    int _cols;
    int _seed;
    double _probability;
    std::vector<std::vector<char>> _world;
    std::vector<std::vector<char>> _worldCopy;
    MPI_Comm _cart;
    MPI_Datatype _colType;
    MPI_Datatype _rowType;
    MPI_Datatype _cornerType;
    int _ndims = 2;
    int _dims[2];
    int _period[2];

    int _upperRightRank;
    int _lowerLeftRank;
    int _lowerRightRank;
    int _upperLeftRank;

    char getRandomValue(int row_id, int col_id);
    int countAlive();
    int countNeighbours(int i, int j);
    int countNeighbours_parallel(int i, int j);

    int plus(int x, int m);
    int minus(int x, int m);
};

#endif // GAMEOFLIFE_H
