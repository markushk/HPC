#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <mpi.h>


char get_random_value(const int row_id, const int col_id, const int n,
			 const int seed) {
    char r = '0';
    int my_seed = seed + row_id * n + col_id;
    //printf("my_seed: %d\n", my_seed);
    srand(my_seed);
    //rand();
    //printf("rand=%d\n", rand());
    rand();
    r = (rand() % 2) ? '1' : '0';
    return r;
}

void initial_configuration(std::vector<std::vector<char>> &world, int rows, int cols, int seed) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            world[i][j] = get_random_value(i, j, cols, seed);
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
    std::cout << "number of alive cells: " << alive << std::endl;
    std::cout << "number of dead cells: " << rows*cols - alive << std::endl;
}

void print_field(std::vector<std::vector<char>> &world, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            if (world[i][j] == '1') std::cout << " X";
            else std::cout << " O";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

int plus(int x, int m) {
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
    if (argc != 6) {
        std::cerr << "Usage: " << argv[0] << " <generations> <rows> <cols> <type> <repetitions>\n";
        return EXIT_FAILURE;
    }
    int generations = std::stoi(argv[1]);
    int rows = std::stoi(argv[2]);
    int cols = std::stoi(argv[3]);
    int type = std::stoi(argv[4]);
    int repetitions = std::stoi(argv[5]);

    for (int i=0; i<repetitions; i++) {
        std::vector<std::vector<char>> world(rows, std::vector<char>(cols));
        std::vector<std::vector<char>> world_copy(rows, std::vector<char>(cols));
        initial_configuration(world, rows, cols, type);
        world_copy = world;
        std::cout << "initial configuration: " << std::endl;
        print_field(world, rows, cols);
        print_status(world, rows, cols);
        // print_field(world_copy, rows, cols);

        double start_time = MPI_Wtime();
        life(world, world_copy, generations, rows, cols);
        double end_time = MPI_Wtime();
        double elapsed_time = end_time - start_time;

        std::cout << "\n\nfinal configuration: " << std::endl;
        print_field(world, rows, cols);

        print_status(world, rows, cols);
        std::cout << "It took " << elapsed_time * 1e6 << " microseconds." << std::endl;
    }

    MPI_Finalize();
    return 0;
}