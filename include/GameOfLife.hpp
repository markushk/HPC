#ifndef GAMEOFLIFE_H
#define GAMEOFLIFE_H

#include <vector>
#include <mpi.h>
#include <iostream>




template <class T>
class matrix2D { 
    std::vector<T> data;
    int columns;
public:
    T &operator()(int x, int y) {
        return data[y * columns + x];
    }

    // matrix2D(int x, int y) : data(x*y), columns(x) {}
    matrix2D(int x, int y, const T &initialValue = T()) : data(x * y, initialValue), columns(x) {}

};

class GameOfLife {
public:
    GameOfLife(int rows, int cols, int seed, double probability,int px,int py);
    ~GameOfLife();

    void initialConfiguration();
    void printStatus();
    void printFieldAnimated();
    void printField();
    void printWholeWorld(int all_rows, int all_cols);
    void printWholeWorldGhost(int all_rows, int all_cols);
    //void runLifep2p(int generations);
    //void runLifep2pCorners(int generations);
    void runLife(int generations, std::string method);
    void init_mpi();
    void finalizempi();
    void testCommunication();
    void pointToPoint();
    void gatherMatrix(int all_rows, int all_cols);
    void gatherMatrixGhost(int all_rows, int all_cols);
    void exchangePointToPoint();
    void exchangePointToPointCorners();
    void exchangeCollective();
    void printFieldAll();
    void printStatusAll(int all_rows, int all_cols);
    void backupWorld();
    void restoreWorld();
    //matrix2D<char> _world;
    //matrix2D<char> _worldCopy;
    //matrix2D<char> _wholeWorld;
    //matrix2D<char> _wholeWorldGhost;
    //matrix2D<char> _origWorld;



private:
    int _rows;
    int _cols;
    int _seed;
    double _probability;
    int _px;
    int _py;
     matrix2D<char> _world;
     matrix2D<char> _worldCopy;
     matrix2D<char> _wholeWorld;
     matrix2D<char> _wholeWorldGhost;
     matrix2D<char> _origWorld;
    MPI_Comm _cart;
    MPI_Comm _coll;
    MPI_Datatype _colType;
    MPI_Datatype _rowType;
    MPI_Datatype _fullColType;
    MPI_Datatype _fullRowType;
    MPI_Datatype _cornerType;
    MPI_Datatype _ROW, _COL, _COR;
    int _ndims = 2;
    int _dims[2];
    int _period[2];
    int _coords[2];

    int _upperRightRank;
    int _lowerLeftRank;
    int _lowerRightRank;
    int _upperLeftRank;
    int _top;
    int _bot;
    int _left;
    int _right;
    int _neighbors[8];
    int sendcount[8], recvcount[8];
    MPI_Aint senddisp[8], recvdisp[8];
    MPI_Datatype sendtype[8], recvtype[8];
    MPI_Request _req[4];
    MPI_Status _statuses[4];

    char getRandomValue(int row_id, int col_id);
    int countAlive();
    int countAliveAll(int all_rows, int all_cols);
    int countNeighbours(int i, int j);

    int plus(int x, int m);
    int minus(int x, int m);
    int c(int x, int y);
};

#endif // GAMEOFLIFE_H
