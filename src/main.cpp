#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include "GameOfLife.hpp"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    if (argc != 9) {
      std::cerr << "Usage: " << argv[0] << " <generations> <rows> <cols> <seed> <probability(%)> <repetitions> <px> <py>\n";
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
    int rank, all_cols, all_rows;
    int p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    all_rows=rows;
    all_cols=cols;
    rows=rows/py;
    cols=cols/px;




    for (int i = 0; i < repetitions; i++) {

        GameOfLife game(rows, cols, seed, probability / 100,px,py);
        game.initialConfiguration();
        /*game.gatherMatrix(all_rows, all_cols);
        game.printField();
        if (rank==0) {
            game.printWholeWorld(all_rows, all_cols);
            //game.printStatusAll(all_rows, all_cols);
            //game.countAliveAll(all_rows, all_cols);
        }*/
        //game.printField();

        //std::cout << "Initial configuration: " << std::endl;
        //game.initialConfiguration();
        //game.printFieldAll();
        //game.exchangePointToPoint();
        //game.printFieldAll();
        //game.printFieldAll();
        //game.exchangePointToPoint();



        /*for (int i = 0; i < p; ++i) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == i) {
                std::cout << "Rank " << rank << " init Matrix:" << std::endl;
                game.printField();
            }
        }*/

        //game.printFieldAll();
        // game.exchangePointToPoint();
        // for (int i = 0; i < p; ++i) {
        //     MPI_Barrier(MPI_COMM_WORLD);
        //     if (rank == i) {
        //         std::cout << "Rank " << rank << " Exchanged Matrix:" << std::endl;
        //         game.printFieldAll();
        //     }
        // }
        //game.printFieldAll();
        // game.printStatus();

        double start_time = MPI_Wtime();
        game.runLifep2p(generations);
        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;
        //std::cout << "\n\nFinal configuration: " << std::endl;
        //game.printField();
        //game.printStatus();
        //std::cout << "It took " << elapsed_time * 1e6 << " microseconds." << std::endl;
        //game.create_mpi();
        //game.testCommunication();
        //game.gatherMatrix();*/
        //game.runLife(generations);
        /*for (int i = 0; i < p; ++i) {
            MPI_Barrier(MPI_COMM_WORLD);
            if (rank == i) {
                std::cout << "Rank " << rank << " Final configuration:" << std::endl;
                game.printField();
            }
        }*/
        //game.printFieldAll();
        game.gatherMatrix(all_rows, all_cols);
        if (rank==0) {
            game.printWholeWorld(all_rows, all_cols);
            game.printStatusAll(all_rows, all_cols);
            //game.countAliveAll(all_rows, all_cols);
        }
    }

    MPI_Finalize();
    return 0;
}
