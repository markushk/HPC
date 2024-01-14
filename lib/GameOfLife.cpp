#include "GameOfLife.hpp"
#include <iostream>
#include <cstdlib>
#include <thread>
#include <chrono>
#include <mpi.h>
#include <algorithm>
#include <iterator>

GameOfLife::GameOfLife(int rows, int cols, int seed, double probability,int px,int py)
    : _rows(rows + 2), _cols(cols + 2), _seed(seed), _probability(probability),_px(px),_py(py),
      _world(_rows, _cols,'#'),
       _worldCopy(_rows, _cols,'#'),
      _wholeWorld(rows*py, cols*px,'#'),
      _wholeWorldGhost(rows*py+2*py, cols*px+2*px,'#'),
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
            //_world(i,j) = getRandomValue((i-1) + coords[1] * (_rows-2), (j-1) + coords[0] * (_cols-2));
            _world(i,j) = '0'+((i-1) + coords[1] * (_rows-2))*_rows+((j-1) + coords[0] * (_cols-2));
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

int GameOfLife::c(int x, int y) {
    return y * _rows + x;
}

// 0 5 10 15 20 25 30 35
// 1 6 11 16 21 26 31 36
// 2 7 12 17 22 27 32 37
// 3 8 13 18 23 28 33 38
// 4 9 14 19 24 29 34 39


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
            if (_wholeWorld(i,j) == '1')
                alive++;
        }
    }
    return alive;
}

void GameOfLife::printFieldAnimated() {
    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {
            if (_world(i,j) == '1')
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
            if (_world(i,j) == '1')
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
            if (_wholeWorld(i,j) == '1')
                std::cout << " X";
            else
                std::cout << " .";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void GameOfLife::printWholeWorldGhost(int all_rows, int all_cols) {
    all_rows=all_rows+2*_py;
    all_cols=all_cols+2*_px;
    int counter = 0;
    for (int i = 0; i < all_rows; ++i) {
        counter=0;
        for (int j = 0; j < all_cols; ++j) {
            /*if (j % _cols == 0 && j != 0 ) {
                std::cout << "  ";

            }
            if (i % _rows == 0 && i != 0 && counter < 1 ) {
                std::cout << "\n";
                counter++;
                //continue;
            }
            if ((i % _rows != 0 || j % _cols != 0) || i == 0 || j == 0) {

                if (_wholeWorldGhost(i,j) == '1')
                        std::cout << " X";
                    else
                        std::cout << " .";
            }*/
            if (j % _cols == 0 && j != 0 ) {
                std::cout << "  ";

            }
            if (i % _rows == 0 && i != 0 && counter < 1 ) {
                std::cout << "\n";
                counter++;
                //continue;
            }
            if ((i % _rows != 0 || j % _cols != 0) || i == 0 || j == 0) {

                std::cout << _wholeWorldGhost(i,j) << " ";
            }
            
                
                
        }
        std::cout << std::endl;
    }
    //counter=0;
    std::cout << std::endl << std::endl ;
}

/*void GameOfLife::printWholeWorldGhost(int all_rows, int all_cols) {
    all_rows=all_rows+2*_py;
    all_cols=all_cols+2*_px;
    for (int i = 0; i < all_rows; ++i) {
        for (int j = 0; j < all_cols; ++j) {
            if (_wholeWorldGhost(i,j) == '1')
                std::cout << " X";
            else
                std::cout << " .";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}*/


void GameOfLife::printFieldAll() {
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            std::cout <<" " << _world(i,j);
            // if (_world(i,j) == '1')
            //     std::cout << " X";
            // else
            //     std::cout << " .";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void GameOfLife::runLife(int generations, std::string method="p2p") {
    for (int gen = 0; gen < generations; gen++) {
        if (method == "p2p") {exchangePointToPoint();};
        if (method == "p2pc") {exchangePointToPointCorners();};
        if (method == "coll") {exchangeCollective();};
        for (int i = 1; i < _rows-1; ++i) {
            for (int j = 1; j < _cols-1; ++j) {
                int neighbours = countNeighbours(i, j);

                if (_world(i,j) == '0') {
                    if (neighbours == 3) {
                        _worldCopy(i,j) = '1';
                    } else {
                        _worldCopy(i,j) = '0';
                    }
                } else {
                    if (neighbours != 2 && neighbours != 3) {
                        _worldCopy(i,j) = '0';
                    } else {
                        _worldCopy(i,j) = '1';
                    }
                }
            }
        }
        // std::copy(&_world(0,0), &_world(0,0)+_rows*_cols,&_worldCopy(0,0));
        _world = _worldCopy;
    }
}





int GameOfLife::countAlive() {
    int alive = 0;
    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {
            if (_world(i,j) == '1')
                alive++;
        }
    }
    return alive;
}

int GameOfLife::plus(int x, int m) {
    return (x + 1) % (m);
}

int GameOfLife::minus(int x, int m) {
    return (x - 1 + m) % (m);
}


int GameOfLife::countNeighbours(int i, int j) {
    int number_of_neighbours = 0;
    if (_world(i,j+1) == '1') number_of_neighbours++;     // right
    if (_world(i+1,j+1) == '1') number_of_neighbours++;  // down right
    if (_world(i+1,j) == '1') number_of_neighbours++;    // down
    if (_world(i+1,j-1) == '1') number_of_neighbours++;  // down left
    if (_world(i,j-1) == '1') number_of_neighbours++;    // left
    if (_world(i-1,j-1) == '1') number_of_neighbours++;   // up left
    if (_world(i-1,j) == '1') number_of_neighbours++;     // up
    if (_world(i-1,j+1) == '1') number_of_neighbours++;  // up right
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
    dims[0]=_px;
    dims[1]=_py;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, period, 0, &_cart);
    MPI_Comm_rank(_cart, &rank);

    MPI_Type_vector(_cols, 1, _rows, MPI_CHAR, &_fullRowType);
    MPI_Type_commit(&_fullRowType);
    MPI_Type_contiguous(_rows, MPI_CHAR, &_fullColType);
    MPI_Type_commit(&_fullColType);

    MPI_Type_vector(_cols-2, 1, _rows, MPI_CHAR, &_rowType);
    MPI_Type_commit(&_rowType);
    MPI_Type_contiguous(_rows-2, MPI_CHAR, &_colType);
    MPI_Type_commit(&_colType);


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
    //MPI_Cart_shift(_cart, 1, 1, &_top, &_bot);
    //MPI_Cart_shift(_cart, 0, 1, &_left, &_right);
    MPI_Cart_shift(_cart, 1, 1, &_top, &_bot);
    MPI_Cart_shift(_cart, 0, 1, &_left, &_right);

    MPI_Comm_rank(_cart, &rank);
    MPI_Cart_coords(_cart, rank, 2, _coords);

    // 0 1 2
    // 3 4 5
    // 6 7 8

    // 0 3 6
    // 1 4 7
    // 2 5 8

    // rank 4:
    // 0 5 10 15 20
    // 1 6 11 16 21
    // 2 7 12 17 22
    // 3 8 13 18 23
    // 4 9 14 19 24

    // 0 1 2 3 4
    // 5 6 7 8 9
    // 10 11 12 13 14


    _neighbors[0] = _upperLeftRank;
    _neighbors[1] = _upperRightRank;
    _neighbors[2] = _lowerLeftRank;
    _neighbors[3] = _lowerRightRank;
    _neighbors[4] = _left;
    _neighbors[5] = _right;
    _neighbors[6] = _top;
    _neighbors[7] = _bot;

    MPI_Dist_graph_create_adjacent(_cart, 8, _neighbors, MPI_UNWEIGHTED, 
                                          8, _neighbors, MPI_UNWEIGHTED, 
                                          MPI_INFO_NULL, 0, &_coll);
    

    // MPI_Type_contiguous(_cols-2,MPI_CHAR,&_ROW);
    // MPI_Type_commit(&_ROW);
    // MPI_Type_vector(_rows-2,1,_cols,MPI_CHAR,&_COL);
    // MPI_Type_commit(&_COL);
    // MPI_Type_contiguous(1,MPI_CHAR,&_COR);
    // MPI_Type_commit(&_COR);
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
    matrix2D<char> wholeWorld(all_rows,all_cols);
    std::vector<char> partMatrix((_rows-2)*(_cols-2));
    std::vector<char> recv_buf((_rows-2)*(_cols-2));
    for (int i = 1; i < _rows-1; ++i) {
        for (int j = 1; j < _cols-1; ++j) {
            partMatrix[(i-1) * (_cols-2) + j-1]=_world(i,j);
            wholeWorld(i-1,j-1)=_world(i,j);
        }
    }

    if (rank != 0) {
        MPI_Send(&partMatrix[0],(_rows-2)*(_cols-2) , MPI_CHAR, 0, 1, _cart);
    }

    if (rank==0) {
        for (int i = 1; i < p; ++i) {
            MPI_Recv(&recv_buf[0],(_rows-2)*(_cols-2) , MPI_CHAR, i, 1, _cart, MPI_STATUS_IGNORE);
            MPI_Cart_coords(_cart, i, 2, coords);


            for (int i = 0; i < _rows-2; ++i) {
                for (int j = 0; j < _cols-2; ++j) {
                    //std::cout << recv_buf[i * (_cols-2) + j] << " ";
                    wholeWorld(i+coords[1]*(_rows-2),j+coords[0]*(_cols-2))=recv_buf[i * (_cols-2) + j];
                }
            }
        }
        for (int i = 0; i < all_rows; ++i) {
            for (int j = 0; j < all_cols; ++j) {
                _wholeWorld(i,j)=wholeWorld(i,j);
            }
        }
    }
}

void GameOfLife::gatherMatrixGhost(int all_rows, int all_cols) {
    int rank, p;
    MPI_Comm_rank(_cart, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int coords[2];
    matrix2D<char> wholeWorldGhost(all_rows+_py*2,all_cols+_px*2);
    std::vector<char> partMatrix((_rows)*(_cols));
    std::vector<char> recv_buf((_rows)*(_cols));
    for (int i = 0; i < _rows; ++i) {
        for (int j = 0; j < _cols; ++j) {
            partMatrix[(i) * (_cols) + j]=_world(i,j);
            wholeWorldGhost(i,j)=_world(i,j);
        }
    }

    if (rank != 0) {
        MPI_Send(&partMatrix[0],(_rows)*(_cols) , MPI_CHAR, 0, 1, _cart);
    }

    if (rank==0) {
        for (int i = 1; i < p; ++i) {
            MPI_Recv(&recv_buf[0],(_rows)*(_cols) , MPI_CHAR, i, 1, _cart, MPI_STATUS_IGNORE);
            MPI_Cart_coords(_cart, i, 2, coords);


            for (int i = 0; i < _rows; ++i) {
                for (int j = 0; j < _cols; ++j) {
                    //std::cout << recv_buf[i * (_cols-2) + j] << " ";
                    wholeWorldGhost(i+coords[1]*(_rows),j+coords[0]*(_cols))=recv_buf[i * (_cols) + j];
                }
            }
        }
        for (int i = 0; i < all_rows + _py*2; ++i) {
            for (int j = 0; j < all_cols + _px*2; ++j) {
                _wholeWorldGhost(i,j)=wholeWorldGhost(i,j);
            }
        }
    }
}


void GameOfLife::exchangePointToPoint() {

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
    MPI_Irecv(&_world(0,_cols-1), 1, _fullColType, _right, 7, _cart, &req[0]);
    MPI_Irecv(&_world(0,0), 1, _fullColType, _left, 8, _cart, &req[1]);
    MPI_Isend(&_world(0,1), 1, _fullColType, _left, 7, _cart, &req[3]);
    MPI_Isend(&_world(0,_cols-2), 1, _fullColType, _right, 8, _cart, &req[2]);
    MPI_Waitall(2, req, statuses);

    MPI_Irecv(&_world(_rows-1,0), 1, _fullRowType, _bot, 5, _cart, &req[1]);
    MPI_Irecv(&_world(0,0), 1, _fullRowType, _top, 6, _cart, &req[0]);
    MPI_Isend(&_world(1,0), 1, _fullRowType, _top, 5, _cart, &req[2]);
    MPI_Isend(&_world(_rows-2,0), 1, _fullRowType, _bot, 6, _cart, &req[3]);
    MPI_Waitall(2, req, statuses);

    // Send to upper left
    /*MPI_Isend(&_world(1,1), 1, MPI_CHAR, _upperLeftRank, 1, _cart, &send_request[0]);
    MPI_Irecv(&_world(0,0), 1, MPI_CHAR, _upperLeftRank, 0, _cart, &recv_request[0]);

    // Send to upper right
    MPI_Isend(&_world(1,_cols-2), 1, MPI_CHAR, _upperRightRank, 2, _cart, &send_request[1]);
    MPI_Irecv(&_world(0,_cols-1), 1, MPI_CHAR, _upperRightRank, 4, _cart, &recv_request[1]);

    // Send to lower left
    MPI_Isend(&_world(_rows-2,1), 1, MPI_CHAR, _lowerLeftRank, 4, _cart, &send_request[2]);
    MPI_Irecv(&_world(_rows-1,0), 1, MPI_CHAR, _lowerLeftRank, 2, _cart, &recv_request[2]);

    // Send to lower right
    MPI_Isend(&_world(_rows-2,_cols-2), 1, MPI_CHAR, _lowerRightRank, 0, _cart, &send_request[3]);
    MPI_Irecv(&_world(_rows-1,_cols-1), 1, MPI_CHAR, _lowerRightRank, 1, _cart, &recv_request[3]);*/

    // Wait for all operations to complete
    //MPI_Waitall(4, send_request, send_status);
    //MPI_Waitall(4, recv_request, recv_status);
    //MPI_Waitall(4, req, statuses);

}

void GameOfLife::exchangePointToPointCorners() {

    //std::cout << "\n";
     
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


    MPI_Request send_request[4], recv_request[4];
    MPI_Status send_status[4], recv_status[4];
    MPI_Request req[4];
    MPI_Status statuses[4];

    MPI_Type_vector(_cols-2, 1, _rows, MPI_CHAR, &_rowType);
    MPI_Type_commit(&_rowType);
    MPI_Type_contiguous(_rows-2, MPI_CHAR, &_colType);
    MPI_Type_commit(&_colType);
    MPI_Type_contiguous(1, MPI_CHAR, &_cornerType);
    MPI_Type_commit(&_cornerType);

    MPI_Irecv(&_world(1,_cols-1), 1, _colType, _right, 7, _cart, &req[0]);
    MPI_Irecv(&_world(1,0), 1, _colType, _left, 8, _cart, &req[1]);
    MPI_Isend(&_world(1,1), 1, _colType, _left, 7, _cart, &req[3]);
    MPI_Isend(&_world(1,_cols-2), 1, _colType, _right, 8, _cart, &req[2]);
    

    MPI_Irecv(&_world(_rows-1,1), 1, _rowType, _bot, 5, _cart, &req[1]);
    MPI_Irecv(&_world(0,1), 1, _rowType, _top, 6, _cart, &req[0]);
    MPI_Isend(&_world(1,1), 1, _rowType, _top, 5, _cart, &req[2]);
    MPI_Isend(&_world(_rows-2,1), 1, _rowType, _bot, 6, _cart, &req[3]);
    

    MPI_Waitall(4, req, statuses);
    // Send to upper left
    MPI_Isend(&_world(1,1), 1, MPI_CHAR, _upperLeftRank, 1, _cart, &send_request[0]);
    MPI_Irecv(&_world(0,0), 1, MPI_CHAR, _upperLeftRank, 0, _cart, &recv_request[0]);

    // Send to upper right
    MPI_Isend(&_world(1,_cols-2), 1, MPI_CHAR, _upperRightRank, 2, _cart, &send_request[1]);
    MPI_Irecv(&_world(0,_cols-1), 1, MPI_CHAR, _upperRightRank, 4, _cart, &recv_request[1]);

    // Send to lower left
    MPI_Isend(&_world(_rows-2,1), 1, MPI_CHAR, _lowerLeftRank, 4, _cart, &send_request[2]);
    MPI_Irecv(&_world(_rows-1,0), 1, MPI_CHAR, _lowerLeftRank, 2, _cart, &recv_request[2]);

    // Send to lower right
    MPI_Isend(&_world(_rows-2,_cols-2), 1, MPI_CHAR, _lowerRightRank, 0, _cart, &send_request[3]);
    MPI_Irecv(&_world(_rows-1,_cols-1), 1, MPI_CHAR, _lowerRightRank, 1, _cart, &recv_request[3]);
    // Wait for all operations to complete
    MPI_Waitall(4, send_request, send_status);
    MPI_Waitall(4, recv_request, recv_status);
    MPI_Waitall(4, req, statuses);



    //printFieldAll();

}

void GameOfLife::exchangeCollective() {
    
    
    int t = 8;
    int rank;
    MPI_Comm_rank(_coll, &rank);
    // for (int i = 0; i < 9; ++i) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (rank == i) {
    //         std::cout << "rank " << rank << " before sending:\n";
    //         printFieldAll();
    //     }
    // }

    int sendcount[t], recvcount[t];
    MPI_Aint senddisp[t], recvdisp[t];
    MPI_Datatype sendtype[t], recvtype[t];

    for (int i=0; i<t; i++) {
        sendcount[i] = 1;
        recvcount[i] = 1;
    }

    // order: top left, top right, bottom left, bottom right, left, right, top, bottom,  
    senddisp[0] = c(1, 1);               sendtype[0] =  _cornerType; // top left
    senddisp[1] = c(1, _cols-2);         sendtype[1] =  _cornerType; // top right
    senddisp[2] = c(_rows-2, 1);         sendtype[2] =  _cornerType; // bottom left
    senddisp[3] = c(_rows-2, _cols-2);   sendtype[3] =  _cornerType; // bottom left
    senddisp[4] = senddisp[0];           sendtype[4] =  _colType; // left
    senddisp[5] = senddisp[1];           sendtype[5] =  _colType; // right
    senddisp[6] = senddisp[0];           sendtype[6] =  _rowType; // top
    senddisp[7] = senddisp[2];           sendtype[7] =  _rowType; // bottom
    // senddisp[4] = senddisp[0];           sendtype[4] =  _colType; // left
    // senddisp[5] = senddisp[1];           sendtype[5] =  _colType; // right
    // senddisp[6] = senddisp[0];           sendtype[6] =  _rowType; // top
    // senddisp[7] = senddisp[2];           sendtype[7] =  _rowType; // bottom



    // ghost layers:
    // order: bottom right, bottom left, top right, top left, right, left, bottom, top,
    // recvdisp[0] = c(_rows-1, _cols-1);   recvtype[0] =  _cornerType; // bottom right
    // recvdisp[1] = c(_rows-1, 0);         recvtype[1] =  _cornerType; // bottom left
    // recvdisp[2] = c(0, _cols-1);         recvtype[2] =  _cornerType; // top right
    // recvdisp[3] = c(0, 0);               recvtype[3] =  _cornerType; // top left
    // recvdisp[4] = c(1, _cols-1);         recvtype[4] =  _colType; // right
    // recvdisp[5] = c(1, 0);               recvtype[5] =  _colType; // left
    // recvdisp[6] = c(_rows-1, 1);         recvtype[6] =  _rowType; // bottom
    // recvdisp[7] = c(0, 1);               recvtype[7] =  _rowType; // top
    bool own_neighbor = false;
    for (int i=0;i<t;i++) {
        if (rank==_neighbors[i])
            own_neighbor = true;
    }
    //own_neighbor = true;
    if (own_neighbor==false) {
             recvdisp[0] = c(0, 0);               recvtype[0] =  _cornerType; // top left
             recvdisp[1] = c(0, _cols-1);         recvtype[1] =  _cornerType; // top right
             recvdisp[2] = c(_rows-1, 0);         recvtype[2] =  _cornerType; // bottom left
             recvdisp[3] = c(_rows-1, _cols-1);   recvtype[3] =  _cornerType; // bottom right
             recvdisp[4] = c(1, 0);               recvtype[4] =  _colType; // left
             recvdisp[5] = c(1, _cols-1);         recvtype[5] =  _colType; // right
             recvdisp[6] = c(0, 1);               recvtype[6] =  _rowType; // top
             recvdisp[7] = c(_rows-1, 1);         recvtype[7] =  _rowType; // bottom
    }
    else {
        recvdisp[0] = c(_rows-1, _cols-1);   recvtype[0] =  _cornerType; // bottom right
        recvdisp[1] = c(_rows-1, 0);         recvtype[1] =  _cornerType; // bottom left
        recvdisp[2] = c(0, _cols-1);         recvtype[2] =  _cornerType; // top right
        recvdisp[3] = c(0, 0);               recvtype[3] =  _cornerType; // top left
        recvdisp[4] = c(1, _cols-1);         recvtype[4] =  _colType; // right
        recvdisp[5] = c(1, 0);               recvtype[5] =  _colType; // left
        recvdisp[6] = c(_rows-1, 1);         recvtype[6] =  _rowType; // bottom
        recvdisp[7] = c(0, 1);               recvtype[7] =  _rowType; // top
    }
    

        
     /*if (rank == 2) {
         std::cout << "neighbors: ";
         for (int i=0; i<t; i++) {
             std::cout <<" " << _neighbors[i];
         }
         std::cout << std::endl;
         std::cout << "senddisp: ";
         for (int i=0; i<t; i++) {
             std::cout <<" " << senddisp[i];
         }
         std::cout << std::endl;
         std::cout << "recvdisp: ";
         for (int i=0; i<t; i++) {
             std::cout <<" " << recvdisp[i];
         }
         std::cout << std::endl;
     }
    // byte offsets
    for (int i=0; i<t; i++) {
        senddisp[i] *= sizeof(char);
        recvdisp[i] *= sizeof(char);
    }*/

    MPI_Neighbor_alltoallw(&_world(0,0),sendcount,senddisp,sendtype,
			               &_world(0,0),recvcount,recvdisp,recvtype,_coll);



    // for (int i = 0; i < 9; ++i) {
    //     MPI_Barrier(MPI_COMM_WORLD);
    //     if (rank == i) {
    //         std::cout << "rank " << rank << " after sending:\n";
    //         printFieldAll();
    //     }
    // }
    // std::cout << "my_rank: " << std::endl;
    // printFieldAll();
}
