#include "Game.hpp"
#include <mpi.h>

Game::Game (const int n_proc, const int n_rows, const int n_cols) : n_proc(n_proc), n_rows(n_rows), n_cols(n_cols)  {

    MPI_Dims_create(n_proc, 2, cartDims);
    n_local_rows = n_rows / cartDims[0];
    n_local_cols = n_cols / cartDims[1];
    Grid.resize(localRows, std::vector<int>(localCols, 0));
    nextGrid.resize(localRows, std::vector<int>(localCols, 0));

}

Game::init_grid() {

    
}