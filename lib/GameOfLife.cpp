#include "GameOfLife.hpp"
#include <iostream>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <mpi.h>

GameOfLife::GameOfLife(int rows, int cols, int seed, double probability,int px,int py)
    : _rows(rows + 2), _cols(cols + 2), _seed(seed), _probability(probability),_px(px),_py(py),
      _world(rows + 2, std::vector<char>(cols + 2)), _worldCopy(rows + 2, std::vector<char>(cols + 2)),
      _wholeWorld(rows*py, std::vector<char>(cols*px)),
      _cart(MPI_COMM_WORLD), _colType(MPI_DATATYPE_NULL), _rowType(MPI_DATATYPE_NULL), _cornerType(MPI_DATATYPE_NULL),
      _upperRightRank(-1), _lowerLeftRank(-1), _lowerRightRank(-1), _upperLeftRank(-1) {
    init_mpi();
}

GameOfLife::~GameOfLife() {}

void GameOfLife::initialConfiguration() {
    int coords[2];
    int rank;
    MPI_Comm_rank(_cart, &rank);
    MPI_Cart_coords(_cart, rank, 2, coords);

    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {

            //_world[i][j] = getRandomValue(i + coords[0] * (_rows-2), j + coords[1] * (_cols -2));
            _world[i][j] = getRandomValue((i-1) + coords[1] * (_rows-2), (j-1) + coords[0] * (_cols-2));
            //_world[i][j] = getRandomValue(i, j);
        }
    }
}

char GameOfLife::getRandomValue(int row_id, int col_id) {
    char r = '0';
    int my_seed = _seed + row_id * (_cols - 2)*_px + col_id;
    srand(my_seed);

    double rand_value = static_cast<double>(rand()) / RAND_MAX;

    if (rand_value < _probability) {
        r = '1';
    }

    return r;
}





void GameOfLife::printStatus() {
    int alive = countAlive();
    std::cout << "Number of alive cells: " << alive << std::endl;
    std::cout << "Number of dead cells: " << (_rows - 2) * (_cols - 2) - alive << std::endl;
}

void GameOfLife::printStatusAll(int all_rows, int all_cols) {
    int alive = countAliveAll(all_rows,all_cols);
    std::cout << "Number of alive cells: " << alive << std::endl;
    std::cout << "Number of dead cells: " << (all_rows) * (all_cols) - alive << std::endl;

}

int GameOfLife::countAliveAll(int all_rows, int all_cols) {
    int alive = 0;
    for (int i = 0; i < all_rows; ++i) {
        for (int j = 0; j < all_cols; ++j) {
            if (_wholeWorld[i][j] == '1')
                alive++;
        }
    }
    return alive;
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

void GameOfLife::printWholeWorld(int all_rows, int all_cols) {
    for (int i = 0; i < all_rows; ++i) {
        for (int j = 0; j < all_cols; ++j) {
            if (_wholeWorld[i][j] == '1')
                std::cout << " X";
            else
                std::cout << " .";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void GameOfLife::printFieldAll() {
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
        exchangePointToPoint();
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
    if (_world[i][j+1] == '1') number_of_neighbours++;     // right
    if (_world[i+1][j+1] == '1') number_of_neighbours++;  // down right
    if (_world[i+1][j] == '1') number_of_neighbours++;    // down
    if (_world[i+1][j-1] == '1') number_of_neighbours++;  // down left
    if (_world[i][j-1] == '1') number_of_neighbours++;    // left
    if (_world[i-1][j-1] == '1') number_of_neighbours++;   // up left
    if (_world[i-1][j] == '1') number_of_neighbours++;     // up
    if (_world[i-1][j+1] == '1') number_of_neighbours++;  // up right
    //std::cout << number_of_neighbours;
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
    //MPI_Dims_create(p, ndims, dims);
    //MPI_Dims_create(p, ndims, _dims);
    dims[0]=_px;
    dims[1]=_py;
    //MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, 0, &_cart);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, 0, &_cart);
    MPI_Comm_rank(_cart, &rank);
    
    MPI_Type_vector(_rows-2, 1, _cols, MPI_CHAR, &_colType);
    MPI_Type_commit(&_colType);
    MPI_Type_contiguous(_cols - 2, MPI_CHAR, &_rowType);
    MPI_Type_commit(&_rowType);
    MPI_Type_contiguous(1, MPI_CHAR, &_cornerType);
    MPI_Type_commit(&_cornerType);

    int coords[2];
    MPI_Cart_coords(_cart, rank, 2, coords);
    
    int upperRightCoords[2] = {(coords[0] + 1 + dims[0]) % dims[0], (dims[1] + coords[1] - 1) % dims[1]};
    int lowerLeftCoords[2] = {(dims[0] + coords[0] - 1) % dims[0], (coords[1] + 1 + dims[1]) % dims[1]};
    int lowerRightCoords[2] = {(coords[0] + 1) % dims[0], (dims[1] + coords[1] + 1) % dims[1]};
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


void GameOfLife::gatherMatrix(int all_rows, int all_cols) {
    int rank, p;
    MPI_Comm_rank(_cart, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int coords[2];
    std::vector<std::vector<char>> wholeWorld(all_rows, std::vector<char>(all_cols));
    std::vector<char> partMatrix((_rows-2)*(_cols-2));
    std::vector<char> recv_buf((_rows-2)*(_cols-2));
    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {
            partMatrix[(i-1) * (_cols-2) + j-1]=_world[i][j];
            wholeWorld[i-1][j-1]=_world[i][j];
            //std::cout << partMatrix[i * _cols + j] << "\n";
        }
    }
    if (rank != 0) {
        MPI_Send(&partMatrix[0],(_rows-2)*(_cols-2) , MPI_CHAR, 0, 1, _cart);
    }
    //MPI_Barrier(_cart);
    //std::cout << "rank: " << p << "\n";
    //std::cout << "all cols: " << all_cols << "\n";
    //std::cout << "all rows: " << all_rows << "\n";


    //MPI_Barrier(_cart);
    if (rank==0) {
        for (int i = 1; i < p; ++i) {
            //std::cout << "getting rank: "<< i << "\n";
            MPI_Recv(&recv_buf[0],(_rows-2)*(_cols-2) , MPI_CHAR, i, 1, _cart, MPI_STATUS_IGNORE);
            MPI_Cart_coords(_cart, i, 2, coords);
            //std::cout << "coords[0]: " << coords[0] << "\n";
            //std::cout << "coords[1]: " << coords[1] << "\n";

            for (int i = 0; i < _rows-2; ++i) {
                for (int j = 0; j < _cols-2; ++j) {
                    //std::cout << recv_buf[i * (_cols-2) + j] << " ";
                    wholeWorld[i+coords[1]*(_rows-2)][j+coords[0]*(_cols-2)]=recv_buf[i * (_cols-2) + j];
                }
                //std::cout << "\n";
            }
        }
        //_wholeWorld=wholeWorld;
        //std::cout << recv_buf[(_cols-2)*(_rows-2)-1] << " ";
        for (int i = 0; i < all_rows; ++i) {
            for (int j = 0; j < all_cols; ++j) {
                _wholeWorld[i][j]=wholeWorld[i][j];
                //std::cout << _wholeWorld[i][j] << " ";
                //_wholeWorld[i][j]=wholeWorld[i][j];
            }
            //std::cout << "\n";

        }


    }


}

void GameOfLife::exchangePointToPoint() {
    int top, bot, left, right, rank;
    int coords[2];

    std::vector<char> sendColLeft(_rows-2);
    std::vector<char> sendColRight(_rows-2);

    

    for (int i = 1; i < _rows-1; ++i) {
        sendColLeft[i-1] = _world[i][1];
        sendColRight[i-1] = _world[i][_cols-2];
    }
        for (int i = 0; i < _rows-2; ++i) {
        //sendColLeft[i-1] = _world[i][1];
        //sendColRight[i-1] = _world[i][_cols-2];
        //std::cout << _world[i][1] << "\n";
        //std::cout << test_send[i][1] << "\n";
        //std::cout << sendColLeft[i] << "\n";
        //_world[_rows][i] = recvRowBot[i];
    }
    //std::cout << "\n";
    
    std::vector<char> rec





    //std::vector<std::vector<char>> world(all_rows, std::vector<char>(all_cols));




    t, rank, 2, coords);

    
    /*if (rank == 9) {
        std::cout << "coord x: " << coords[0] << "\n";
        std::cout << "coord y: " << coords[1] << "\n";
        std::cout << "top: " << top << "\n";
        std::cout << "bot: " << bot << "\n";
        std::cout << "left: " << left << "\n";
        std::cout << "right: " << right << "\n";
        std::cout << "topright: " << _upperRightRank << "\n";
        std::cout << "topleft: " << _upperLeftRank << "\n";
        std::cout << "botleft: " << _lowerLeftRank << "\n";
        std::cout << "botright: " << _lowerRightRank << "\n";

    }*/
    MPI_Request req[4];
    MPI_Status statuses[4];
    if (_py < -3)
    {
        MPI_Irecv(&_world[0][1], 1, _rowType, bot, 1, _cart, &req[1]);
        MPI_Irecv(&_world[_rows-1][1], 1, _rowType, top, 0, _cart, &req[0]);
        MPI_Isend(&_world[1][1], 1, _rowType, top, 0, _cart, &req[2]);
        MPI_Isend(&_world[_rows-2][1], 1, _rowType, bot, 1, _cart, &req[3]);
    }
    else {
        MPI_Irecv(&_world[_rows-1][1], 1, _rowType, bot, 0, _cart, &req[1]);
        MPI_Irecv(&_world[0][1], 1, _rowType, top, 1, _cart, &req[0]);
        MPI_Isend(&_world[1][1], 1, _rowType, top, 0, _cart, &req[2]);
        MPI_Isend(&_world[_rows-2][1], 1, _rowType, bot, 1, _cart, &req[3]);
    }

    MPI_Waitall(4, req, statuses);

    //MPI_Wait(&send_request[0], &send_status[0]);
    //MPI_Wait(&recv_request[0], &recv_status[0]);
    //MPI_Wait(&send_request[1], &send_status[1]);
    //MPI_Wait(&recv_request[1], &recv_status[1]);
    //MPI_Waitall(2, send_request, send_status);
    //MPI_Waitall(2, recv_request, recv_status);


    MPI_Irecv(&recvColLeft[0], _rows-2, MPI_CHAR, right, 0, _cart, &req[0]);
    MPI_Irecv(&recvColRight[0], _rows-2, MPI_CHAR, left, 1, _cart, &req[1]);
    MPI_Isend(&sendColLeft[0], _rows-2, MPI_CHAR, left, 0, _cart, &req[3]);
    MPI_Isend(&sendColRight[0], _rows-2, MPI_CHAR, right, 1, _cart, &req[2]);
    //MPI_Isend(&sendColLeft[0], _rows-2, MPI_CHAR, left, 1, _cart, &req[3]);
    /*MPI_Irecv(&recvColLeft[0], _rows-2, MPI_CHAR, left, 2, _cart, &req[0]);
    MPI_Irecv(&recvColRight[0], _rows-2, MPI_CHAR, right, 3, _cart, &req[1]);
    MPI_Isend(&sendColRight[0], _rows-2, MPI_CHAR, right, 3, _cart, &req[2]);
    MPI_Isend(&sendColLeft[0], _rows-2, MPI_CHAR, left, 2, _cart, &req[3]);*/



    MPI_Waitall(4, req, statuses);
    // Send to upper left
    MPI_Isend(&_world[1][1], 1, MPI_CHAR, _upperLeftRank, 1, _cart, &send_request[0]);
    MPI_Irecv(&_world[0][0], 1, MPI_CHAR, _upperLeftRank, 0, _cart, &recv_request[0]);

    // Send to upper right
    MPI_Isend(&_world[1][_cols-2], 1, MPI_CHAR, _upperRightRank, 2, _cart, &send_request[1]);
    MPI_Irecv(&_world[0][_cols-1], 1, MPI_CHAR, _upperRightRank, 4, _cart, &recv_request[1]);

    // Send to lower left
    MPI_Isend(&_world[_rows-2][1], 1, MPI_CHAR, _lowerLeftRank, 4, _cart, &send_request[2]);
    MPI_Irecv(&_world[_rows-1][0], 1, MPI_CHAR, _lowerLeftRank, 2, _cart, &recv_request[2]);

    // Send to lower right
    MPI_Isend(&_world[_rows-2][_cols-2], 1, MPI_CHAR, _lowerRightRank, 0, _cart, &send_request[3]);
    MPI_Irecv(&_world[_rows-1][_cols-1], 1, MPI_CHAR, _lowerRightRank, 1, _cart, &recv_request[3]);

    // Wait for all operations to complete
    MPI_Waitall(4, send_request, send_status);
    MPI_Waitall(4, recv_request, recv_status);


    for (int i = 1; i < _rows-1; ++i) {
        _world[i][_cols-1] = recvColLeft[i-1];  
        _world[i][0] = recvColRight[i-1];
    }

}
