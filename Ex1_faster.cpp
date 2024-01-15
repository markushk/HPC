#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <chrono>
#include <thread>


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

char get_random_value(const int row_id, const int col_id, const int n, const int seed, double probability) {
    char r = '0';
    int my_seed = seed + row_id * n + col_id;
    //std::cout << "seed: " << my_seed <<"\n";
    srand(my_seed);

    double rand_value = static_cast<double>(rand()) / RAND_MAX;

    if (rand_value < probability) {
        r = '1';
    }

    return r;
}

void initial_configuration(matrix2D<char> &world, int rows, int cols, int seed, double probability) {
    for (int i = 1; i < rows-1; ++i) {
        for (int j = 1; j < cols-1; ++j) {
            world(i,j) = get_random_value(i-1, j-1, cols-2, seed, probability);
            //std::cout << "i,j: " << i << " " << j << "\n";
        }
    }
}



int count_alive(matrix2D<char> &world, int rows, int cols) {
    int alive = 0;
    for (int i = 1; i < rows-1; ++i) {
        for (int j = 1; j < cols-1; ++j) {
            if (world(i,j) == '1') alive++;
        }
    }
    return alive;
}

void print_status(matrix2D<char> &world, int rows, int cols) {
    int alive = count_alive(world, rows, cols);
    std::cout << "Number of alive cells: " << alive << std::endl;
    std::cout << "Number of dead cells: " << (rows-2)*(cols-2) - alive << std::endl;
}


void print_field(matrix2D<char> &world, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (world(i,j) == '1') std::cout << " X";
            else std::cout << " .";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

// int plus(int x, int m) {
//     //std::cout << (x + 1) % (m);
//     return (x+1)%m;
// }

// int minus(int x, int m) {
//     return (x-1+m)%m;
// }

int count_neighbours(matrix2D<char> &world, int i, int j, int rows, int cols) {
    int number_of_neighbours = 0;
    if (world(i,j+1) == '1') number_of_neighbours++;     // right
    if (world(i+1,j+1) == '1') number_of_neighbours++;  // down right
    if (world(i+1,j) == '1') number_of_neighbours++;    // down
    if (world(i+1,j-1) == '1') number_of_neighbours++;  // down left
    if (world(i,j-1) == '1') number_of_neighbours++;    // left
    if (world(i-1,j-1) == '1') number_of_neighbours++;   // up left
    if (world(i-1,j) == '1') number_of_neighbours++;     // up
    if (world(i-1,j+1) == '1') number_of_neighbours++;  // up right
    return number_of_neighbours;
}

void life(matrix2D<char> &world, matrix2D<char> &world_copy, int generations, int rows, int cols) {
    for (int gen = 0; gen<generations; gen++) {
        // print_field(world, rows, cols);
        world(0,0) = world(rows-2,cols-2);
        world(rows-1,cols-1) = world(1,1);
        world(0,cols-1) = world(rows-2,1);
        world(rows-1,0) = world(1,cols-2);
        for (int i = 1; i < cols-1; ++i) {
            world(0,i) = world(rows-2,i);
            world(rows-1,i) = world(1,i);
        }
        for (int j = 1; j < rows-1; ++j) {
            world(j,0) = world(j,cols-2);
            world(j,cols-1) = world(j,1);
        }
        // print_field(world, rows, cols);
        for (int i = 1; i < rows-1; ++i) {
            for (int j = 1; j < cols-1; ++j) {
                int neighbours = count_neighbours(world, i, j, rows, cols);
                if (world(i,j) == '0') {
                    if (neighbours == 3) {
                        world_copy(i,j) = '1';
                    }
                    else {
                        world_copy(i,j) = '0';
                    }
                }
                else {
                    if (neighbours != 2 && neighbours != 3) {
                        world_copy(i,j) = '0';
                    }
                    else {
                        world_copy(i,j) = '1';
                    }
                }
            }
        }
    world = world_copy;
    }
}



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

    rows+=2;
    cols+=2;

    matrix2D<char> world(rows,cols);
    matrix2D<char> world_copy(rows,cols);
    matrix2D<char> world_orig(rows,cols);
    initial_configuration(world_orig, rows, cols, seed, probability/100);

    for (int i=0; i<repetitions; i++) {
        world = world_orig;
        world_copy = world;
        // std::cout << "\n\ninitial configuration: " << std::endl;
        // print_field(world, rows, cols);

        double start_time = MPI_Wtime();
        life(world, world_copy, generations, rows, cols);
        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;

        // std::cout << "\n\nfinal configuration: " << std::endl;
        // print_field(world, rows, cols);

        print_status(world, rows, cols);
        std::cout << "It took " << elapsed_time << " seconds." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
