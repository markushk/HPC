#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include "GameOfLife.hpp"
#include <algorithm>
#include <math.h>

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    if (argc != 10) {
      std::cerr << "Usage: " << argv[0] << " <generations> <rows> <cols> <seed> <probability(%)> <repetitions> <px> <py> <method>\n";
        return EXIT_FAILURE;
    }
    int generations = std::stoi(argv[1]);
    int rows = std::stoi(argv[2]);
    int cols = std::stoi(argv[3]);
    int seed = std::stoi(argv[4]);
    double probability = std::stod(argv[5]);
    int repetitions = std::stoi(argv[6]);
    int px = std::stoi(argv[7]);
    int py = std::stoi(argv[8]);
    std::string method = std::string(argv[9]);



    int rank, all_cols, all_rows;
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    all_rows=rows;
    all_cols=cols;
    rows=rows/py;
    cols=cols/px;



    std::vector<double> times(repetitions);
    for (int i = 0; i < repetitions; i++) {



        GameOfLife game(rows, cols, seed, probability / 100,px,py);
        game.initialConfiguration();
        int debug = 0;
        if (debug ==1) {
            std::cout << "debug enabled \n";
            game.gatherMatrix(all_rows, all_cols);
            if (rank==0) {
                game.printWholeWorld(all_rows, all_cols);
            }
            game.gatherMatrixGhost(all_rows, all_cols);

            if (rank==0) {
                game.printWholeWorldGhost(all_rows, all_cols);
            }
            game.exchangePointToPoint();

            game.gatherMatrixGhost(all_rows, all_cols);
            if (rank==0) {
                game.printWholeWorldGhost(all_rows, all_cols);
            }
            game.printFieldAll();
            game.gatherMatrix(all_rows, all_cols);
            if (rank==0) {
                game.printWholeWorld(all_rows, all_cols);
                game.printStatusAll(all_rows, all_cols);

            }

            }

        MPI_Barrier(MPI_COMM_WORLD);
        double start_time = MPI_Wtime();
        game.runLife(generations,std::string(method));
        MPI_Barrier(MPI_COMM_WORLD);
        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;
        times[i]=elapsed_time;

    }
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank==0) {
    double kernel_sum = 0.0;
        for (int k = 0; k < repetitions; k++) {
            kernel_sum += times[k];
        }
        double average = kernel_sum / repetitions;
        double error = 0.0;
        for (int k = 0; k < repetitions; k++) {
            error += pow(times[k] - average, 2);
        }
        error = sqrt(error / (repetitions - 1));

        printf("%d %d %d %d %.8f %.8f\n", px, py, all_rows, all_cols, average, error);
    }

    MPI_Finalize();
    return 0;
}
