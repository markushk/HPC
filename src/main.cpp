#include <cstdlib>
#include <iostream>
#include <mpi.h>
#include "GameOfLife.hpp"

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    if (argc != 7) {
        std::cerr << "Usage: " << argv[0] << " <generations> <rows> <cols> <seed> <probability(%)> <repetitions>\n";
        return EXIT_FAILURE;
    }
    int generations = std::stoi(argv[1]);
    int rows = std::stoi(argv[2]);
    int cols = std::stoi(argv[3]);
    int seed = std::stoi(argv[4]);
    double probability = std::stod(argv[5]);
    int repetitions = std::stoi(argv[6]);
    int rank, size;
    int p;



    for (int i = 0; i < repetitions; i++) {

        GameOfLife game(rows, cols, seed, probability / 100);
        /*game.initialConfiguration();

        std::cout << "Initial configuration: " << std::endl;
        game.printField();
        game.printStatus();

        double start_time = MPI_Wtime();
        game.runLife(generations);
        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;
        std::cout << "\n\nFinal configuration: " << std::endl;
        game.printField();
        game.printStatus();
        std::cout << "It took " << elapsed_time * 1e6 << " microseconds." << std::endl;*/
        //game.create_mpi();
        game.testCommunication();
    }

    MPI_Finalize();
    return 0;
}
