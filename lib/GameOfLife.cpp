#include "GameOfLife.hpp"
#include <iostream>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <mpi.h>

GameOfLife::GameOfLife(int rows, int cols, int seed, double probability)
    : _rows(rows + 2), _cols(cols + 2), _seed(seed), _probability(probability),
      _world(rows + 2, std::vector<char>(cols + 2)), _worldCopy(rows + 2, std::vector<char>(cols + 2)),
      _cart(MPI_COMM_WORLD), _colType(MPI_DATATYPE_NULL), _rowType(MPI_DATATYPE_NULL), _cornerType(MPI_DATATYPE_NULL),
      _upperRightRank(-1), _lowerLeftRank(-1), _lowerRightRank(-1), _upperLeftRank(-1) {
    init_mpi();
}

GameOfLife::~GameOfLife() {}

void GameOfLife::initialConfiguration() {
    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {

            _world[i][j] = getRandomValue(i, j);
        }
    }
}

void GameOfLife::printStatus() {
    int alive = countAlive();
    std::cout << "Number of alive cells: " << alive << std::endl;
    std::cout << "Number of dead cells: " << (_rows - 2) * (_cols - 2) - alive << std::endl;
}

void GameOfLife::printFieldAnimated() {
    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {
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
    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {
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
        for (int i = 1; i < _rows-1; ++i) {
            for (int j = 1; j < _cols-1; ++j) {
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
    int my_seed = _seed + (row_id - 1) * (_cols - 2) + col_id -1;
    srand(my_seed);

    double rand_value = static_cast<double>(rand()) / RAND_MAX;

    if (rand_value < _probability) {
        r = '1';
    }

    return r;
}

int GameOfLife::countAlive() {
    int alive = 0;
    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {
            if (_world[i][j] == '1')
                alive++;
        }
    }
    return alive;
}

int GameOfLife::plus(int x, int m) {
    //std::cout << (x + 1) % (m);
    return (x + 1) % (m);
}

int GameOfLife::minus(int x, int m) {
    return (x - 1 + m) % (m);
}



int GameOfLife::countNeighbours(int i, int j) {
    int number_of_neighbours = 0;
    if (_world[i][plus(j, _cols-2)] == '1') number_of_neighbours++;     // right
    if (_world[minus(i, _rows-2)][plus(j, _cols-2)] == '1') number_of_neighbours++;  // down right
    if (_world[minus(i, _rows-2)][j] == '1') number_of_neighbours++;    // down
    if (_world[minus(i, _rows-2)][minus(j, _cols-2)] == '1') number_of_neighbours++;  // down left
    if (_world[i][minus(j, _cols-2)] == '1') number_of_neighbours++;    // left
    if (_world[plus(i, _rows-2)][minus(j, _cols-2)] == '1') number_of_neighbours++;   // up left
    if (_world[plus(i, _rows-2)][j] == '1') number_of_neighbours++;     // up
    if (_world[plus(i, _rows-2)][plus(j, _cols-2)] == '1') number_of_neighbours++;  // up right
    std::cout << number_of_neighbours;
    return number_of_neighbours;
}

int GameOfLife::countNeighbours_parallel(int i, int j) {
    int number_of_neighbours = 0;
    if (_world[i][j+1] == '1') number_of_neighbours++;     // right
    if (_world[i+1][j+1] == '1') number_of_neighbours++;  // down right
    if (_world[i+1][j] == '1') number_of_neighbours++;    // down
    if (_world[i+1][j-1] == '1') number_of_neighbours++;  // down left
    if (_world[i][j-1] == '1') number_of_neighbours++;    // left
    if (_world[i-1][j-1] == '1') number_of_neighbours++;   // up left
    if (_world[i-1][j] == '1') number_of_neighbours++;     // up
    if (_world[i-1][j+1] == '1') number_of_neighbours++;  // up right
    std::cout << number_of_neighbours;
    return number_of_neighbours;
}




void GameOfLife::init_mpi() {
    int p, rank;
    int ndims = 2;
    int dims[ndims];
    int period[ndims];
    for (int i = 0; i < ndims; i++) {
        dims[i] = 0;
        period[i] = 1;
    }
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Dims_create(p, ndims, dims);
    //MPI_Dims_create(p, ndims, _dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, 0, &_cart);
    MPI_Comm_rank(_cart, &rank);
    MPI_Type_vector(_rows - 2, 1, _cols, MPI_CHAR, &_colType);
    MPI_Type_commit(&_colType);
    MPI_Type_contiguous(_cols - 2, MPI_CHAR, &_rowType);
    MPI_Type_commit(&_rowType);
    MPI_Type_contiguous(1, MPI_CHAR, &_cornerType);
    MPI_Type_commit(&_cornerType);

    int coords[2];
    MPI_Cart_coords(_cart, rank, 2, coords);

    int upperRightCoords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] + 1) % dims[1]};
    int lowerLeftCoords[2] = {(coords[0] + 1) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};
    int lowerRightCoords[2] = {(coords[0] + 1) % dims[0], (coords[1] + 1) % dims[1]};
    int upperLeftCoords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};

    MPI_Cart_rank(_cart, upperRightCoords, &_upperRightRank);
    MPI_Cart_rank(_cart, lowerLeftCoords, &_lowerLeftRank);
    MPI_Cart_rank(_cart, lowerRightCoords, &_lowerRightRank);
    MPI_Cart_rank(_cart, upperLeftCoords, &_upperLeftRank);


}


void GameOfLife::finalizempi() {
    if (_colType != MPI_DATATYPE_NULL) {
        MPI_Type_free(&_colType);
    }
    if (_rowType != MPI_DATATYPE_NULL) {
        MPI_Type_free(&_rowType);
    }
    if (_cornerType != MPI_DATATYPE_NULL) {
        MPI_Type_free(&_cornerType);
    }
}


void GameOfLife::pointToPoint() {

}


void GameOfLife::testCommunication() {
    int rank, size;

    int p;
    int ndims = 2;
    int dims[ndims];
    int period[ndims];
    for (int i = 0; i < ndims; i++) {
        dims[i] = 0;
        period[i] = 1;
    }
    MPI_Comm cart;
    MPI_Comm_size(MPI_COMM_WORLD,&p);
    MPI_Dims_create(p, ndims, dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, 0, &cart);

    MPI_Datatype colType, rowType, cornerType;
    MPI_Type_vector(_rows-2, 1, _cols, MPI_CHAR, &colType);
    MPI_Type_commit(&colType);
    MPI_Type_contiguous(_cols-2, MPI_CHAR, &rowType);
    MPI_Type_commit(&rowType);
    MPI_Type_contiguous(1, MPI_CHAR, &cornerType);
    MPI_Type_commit(&cornerType);

    MPI_Comm_rank(cart, &rank);
    MPI_Comm_size(cart, &size);
    int coords[2];
    MPI_Cart_coords(cart, rank, 2, coords);

    //int dims[2];
    //MPI_Cart_get(cart, 2, dims, nullptr);

    int upperRightCoords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] + 1) % dims[1]};
    int lowerLeftCoords[2] = {(coords[0] + 1) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};
    int lowerRightCoords[2] = {(coords[0] + 1) % dims[0], (coords[1] + 1) % dims[1]};
    int upperLeftCoords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};

    int upperRightRank, lowerLeftRank, lowerRightRank, upperLeftRank;

    MPI_Cart_rank(cart, upperRightCoords, &upperRightRank);
    MPI_Cart_rank(cart, lowerLeftCoords, &lowerLeftRank);
    MPI_Cart_rank(cart, lowerRightCoords, &lowerRightRank);
    MPI_Cart_rank(cart, upperLeftCoords, &upperLeftRank);


    std::vector<std::vector<char>> matrix(_rows, std::vector<char>(_cols, rank));


    //std::cout << "Rank " << rank << " Original Matrix:" << std::endl;
    /*for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            matrix[i][j] = rank*j +j;
            std::cout << (int) matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }*/

    int source, dest;
    /*MPI_Cart_shift(cart, 0, 1, &source, &dest);
    MPI_Send(&matrix[0][1], 1, rowType, dest, 0, cart);
    MPI_Recv(&matrix[0][1], 1, rowType, source, 0, cart, MPI_STATUS_IGNORE);

    //MPI_Send(&matrix[_rows-1][1], 1, rowType, dest, 0, cart);
    //MPI_Recv(&matrix[_rows-1][1], 1, rowType, source, 0, cart, MPI_STATUS_IGNORE);
    MPI_Send(&matrix[_rows-1][1], 1, rowType, source, 0, cart);
    MPI_Recv(&matrix[_rows-1][1], 1, rowType, dest, 0, cart, MPI_STATUS_IGNORE);*/

    std::vector<char> sendColLeft(_rows-2);
    std::vector<char> sendColRight(_rows-2);
    for (int i = 1; i < _cols-1; ++i) {
        sendColLeft[i-1] = matrix[i][0];
        sendColRight[i-1] = matrix[i][_cols-1];
        //matrix[_rows][i] = recvRowBot[i];
    }

    std::vector<char> recvColLeft(_rows-2);
    std::vector<char> recvColRight(_rows-2);
    //char recvColLeft[10];
    //char recvColRight[10];

    /*MPI_Cart_shift(cart, 1, 1, &source, &dest);

    MPI_Send(&sendColLeft[0], _rows-2, MPI_CHAR, dest, 0, cart);
    MPI_Recv(&recvColLeft[0], _rows-2, MPI_CHAR, source, 0, cart, MPI_STATUS_IGNORE);
    MPI_Send(&sendColRight[0], _rows-2, MPI_CHAR, source, 0, cart);
    MPI_Recv(&recvColRight[0], _rows-2, MPI_CHAR, dest, 0, cart, MPI_STATUS_IGNORE);*/

// Assuming upperLeftRank, upperRightRank, lowerLeftRank, and lowerRightRank are correctly determined

MPI_Request send_request[4], recv_request[4];
MPI_Status send_status[4], recv_status[4];

MPI_Cart_shift(cart, 0, 1, &source, &dest);
MPI_Isend(&matrix[0][1], 1, rowType, dest, 0, cart, &send_request[0]);
MPI_Irecv(&matrix[0][1], 1, rowType, source, 0, cart, &recv_request[0]);


MPI_Wait(&send_request[0], &send_status[0]);
MPI_Wait(&recv_request[0], &recv_status[0]);

// Send to the right
MPI_Isend(&matrix[_rows-1][1], 1, rowType, source, 0, cart, &send_request[1]);
MPI_Irecv(&matrix[_rows-1][1], 1, rowType, dest, 0, cart, &recv_request[1]);

MPI_Wait(&send_request[1], &send_status[1]);
MPI_Wait(&recv_request[1], &recv_status[1]);

MPI_Cart_shift(cart, 1, 1, &source, &dest);
// Create non-blocking send and receive requests for columns
MPI_Isend(&sendColLeft[0], _rows-2, MPI_CHAR, dest, 0, cart, &send_request[2]);
MPI_Irecv(&recvColLeft[0], _rows-2, MPI_CHAR, source, 0, cart, &recv_request[2]);
MPI_Wait(&send_request[1], &send_status[2]);
MPI_Wait(&recv_request[1], &recv_status[2]);

MPI_Isend(&sendColRight[0], _rows-2, MPI_CHAR, source, 0, cart, &send_request[3]);
MPI_Irecv(&recvColRight[0], _rows-2, MPI_CHAR, dest, 0, cart, &recv_request[3]);
MPI_Wait(&send_request[1], &send_status[3]);
MPI_Wait(&recv_request[1], &recv_status[3]);
// Wait for all non-blocking operations to complete
MPI_Waitall(2, send_request + 2, recv_status + 2);
//MPI_Barrier(cart);

// Send to upper left
MPI_Isend(&matrix[1][1], 1, MPI_CHAR, upperLeftRank, 0, cart, &send_request[0]);
MPI_Irecv(&matrix[0][0], 1, MPI_CHAR, upperLeftRank, 0, cart, &recv_request[0]);

// Send to upper right
MPI_Isend(&matrix[1][_cols-2], 1, MPI_CHAR, upperRightRank, 0, cart, &send_request[1]);
MPI_Irecv(&matrix[0][_cols-1], 1, MPI_CHAR, upperRightRank, 0, cart, &recv_request[1]);

// Send to lower left
MPI_Isend(&matrix[_rows-2][1], 1, MPI_CHAR, lowerLeftRank, 0, cart, &send_request[2]);
MPI_Irecv(&matrix[_rows-1][0], 1, MPI_CHAR, lowerLeftRank, 0, cart, &recv_request[2]);

// Send to lower right
MPI_Isend(&matrix[_rows-2][_cols-2], 1, MPI_CHAR, lowerRightRank, 0, cart, &send_request[3]);
MPI_Irecv(&matrix[_rows-1][_cols-1], 1, MPI_CHAR, lowerRightRank, 0, cart, &recv_request[3]);

// Wait for all operations to complete
MPI_Waitall(4, send_request, send_status);
MPI_Waitall(4, recv_request, recv_status);


    for (int i = 1; i < _cols-1; ++i) {
        //std::cout << (int) recvColLeft[i-1];
        matrix[i][0] = recvColLeft[i-1];  
        matrix[i][_cols-1] = recvColRight[i-1];
        //matrix[_rows][i] = recvRowBot[i];
    }

    // Print in order
    for (int i = 0; i < size; ++i) {
        MPI_Barrier(cart);  // Ensure only one process prints at a time
        if (rank == i) {
            std::cout << "Rank " << rank << " Updated Matrix:" << std::endl;
            for (int i = 0; i < _rows; ++i) {
                for (int j = 0; j < _cols; ++j) {
                    std::cout << (int) matrix[i][j] << " ";
                }
                std::cout << std::endl;
            }
        }
    }


    //MPI_Cart_shift(cart, 1, 1, &source, &dest);
    //MPI_Sendrecv(&matrix[0][0], 1, colType, dest, 0, &recvCol[0], 1, colType, source, 0, cart, MPI_STATUS_IGNORE);
    


    //std::cout << "Rank " << rank << " Received Column:" << std::endl;
    //for (int i = 0; i < _rows-2; ++i) {
    //    std::cout << recvCol[i] << std::endl;
    //}
}
