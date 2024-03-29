#ifndef GAMEOFLIFE_H
#define GAMEOFLIFE_H

#include <vector>
#include <mpi.h>

class GameOfLife {
public:
    GameOfLife(int rows, int cols, int seed, double probability,int px,int py);
    ~GameOfLife();

    void initialConfiguration();
    void printStatus();
    void printFieldAnimated();
    void printField();
    void printWholeWorld(int all_rows, int all_cols);
    void runLife(int generations);
    void init_mpi();
    void finalizempi();
    void testCommunication();
    void pointToPoint();
    void gatherMatrix(int all_rows, int all_cols);
    void exchangePointToPoint();
    void printFieldAll();
    void printStatusAll(int all_rows, int all_cols);



private:
    int _rows;
    int _cols;
    int _seed;
    double _probability;
    int _px;
    int _py;
    std::vector<std::vector<char>> _world;
    std::vector<std::vector<char>> _worldCopy;
    std::vector<std::vector<char>> _wholeWorld;
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
    int countAliveAll(int all_rows, int all_cols);
    int countNeighbours(int i, int j);

    int plus(int x, int m);
    int minus(int x, int m);
};

#endif // GAMEOFLIFE_H
