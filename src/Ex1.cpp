#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <chrono>
#include <thread>
#include <algorithm>
#include <math.h>



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

void initial_configuration(std::vector<std::vector<char>> &world, int rows, int cols, int seed, double probability) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            world[i][j] = get_random_value(i, j, cols, seed, probability);
            //std::cout << "i,j: " << i << " " << j << "\n";
        }
    }
}


int count_alive(std::vector<std::vector<char>> &world, int rows, int cols) {
    int alive = 0;
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (world[i][j] == '1') alive++;
        }
    }
    return alive;
}

void print_status(std::vector<std::vector<char>> &world, int rows, int cols) {
    int alive = count_alive(world, rows, cols);
    std::cout << "Number of alive cells: " << alive << std::endl;
    std::cout << "Number of dead cells: " << rows*cols - alive << std::endl;
}

void print_field_animated(std::vector<std::vector<char>> &world, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (world[i][j] == '1') std::cout << " X";
            else std::cout << " .";
        }
        std::cout << "\n";
    }
    std::cout << "\x1B[" << rows << "A";  // Carriage return to move back to the beginning of the line
    std::this_thread::sleep_for(std::chrono::milliseconds(300));  // Adjust delay as needed
    std::cout.flush();  // Flush the output
}

void print_field(std::vector<std::vector<char>> &world, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (world[i][j] == '1') std::cout << " X";
            else std::cout << " .";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int plus(int x, int m) {
    //std::cout << (x + 1) % (m);
    return (x+1)%m;
}

int minus(int x, int m) {
    return (x-1+m)%m;
}

int count_neighbours(std::vector<std::vector<char>> &world, int i, int j, int rows, int cols) {
    int number_of_neighbours = 0;
    if (world[i][plus(j,cols)] == '1') number_of_neighbours++; // right
    if (world[minus(i,rows)][plus(j,cols)] == '1') number_of_neighbours++; // down right
    if (world[minus(i,rows)][j] == '1') number_of_neighbours++; // down
    if (world[minus(i,rows)][minus(j,cols)] == '1') number_of_neighbours++; // down left
    if (world[i][minus(j,cols)] == '1') number_of_neighbours++; // left
    if (world[plus(i,rows)][minus(j,cols)] == '1') number_of_neighbours++; // up left
    if (world[plus(i,rows)][j] == '1') number_of_neighbours++; // up
    if (world[plus(i,rows)][plus(j,cols)] == '1') number_of_neighbours++; // up right
    //std::cout << number_of_neighbours;
    return number_of_neighbours;
}

void life(std::vector<std::vector<char>> &world, std::vector<std::vector<char>> &world_copy, int generations, int rows, int cols) {
    for (int gen = 0; gen<generations; gen++) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                int neighbours = count_neighbours(world, i, j, rows, cols);
                if (world[i][j] == '0') {
                    if (neighbours == 3) {
                        world_copy[i][j] = '1';
                    }
                }
                else {
                    if (neighbours != 2 && neighbours != 3) {
                        world_copy[i][j] = '0';
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
    std::vector<double> times(repetitions);
    std::vector<std::vector<char>> world(rows, std::vector<char>(cols));
    std::vector<std::vector<char>> world_copy(rows, std::vector<char>(cols));
    std::vector<std::vector<char>> world_orig(rows, std::vector<char>(cols));

    initial_configuration(world_orig, rows, cols, seed, probability/100);
    for (int i=0; i<repetitions; i++) {
        world=world_orig;
        world_copy = world;

        double start_time = MPI_Wtime();
        life(world, world_copy, generations, rows, cols);
        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;
        times[i]=elapsed_time;

    }
        print_field(world, rows, cols);
        print_status(world, rows, cols);
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

        printf(" %d %d %.8f %.8f\n", rows, cols, average, error);


    MPI_Finalize();
    return 0;
}