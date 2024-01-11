#include <cstdlib>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>
#include <mpi.h>
#include <chrono>
#include <thread>

void print_field(std::vector<std::vector<char>> &world, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << " " << world[i][j];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}


char get_random_value(const int row_id, const int col_id, const int n, const int seed, double probability) {
    char r = '0';
    int my_seed = seed + row_id * n + col_id;
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
        }
    }
}


int main(int argc, char *argv[]) {
    //RUN WITH: mpirun --oversubscribe -np 9 ./filename 1 9 9 0 50 1 3 3
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
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);    



    int rank_rows = rows/py + 2;
    int rank_cols = cols / px + 2;

    std::vector<std::vector<char>> world(rank_rows, std::vector<char>(rank_cols,'0'+rank));
    // initial_configuration(world, rank_rows, rank_cols, seed, probability/100);

    if (rank == 4) {
        std::cout << "initial configuration: \n";
        print_field(world, rank_rows, rank_cols);
    }

    // Create a Cartesian communicator for the 2D grid
    MPI_Comm grid_comm;
    int dims[2] = {px, py};
    int periods[2] = {1, 1}; // Periodic grid in both dimensions
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &grid_comm);

    // Determine the coordinates of the current process in the 2D grid
    int coords[2];
    MPI_Cart_coords(grid_comm, rank, 2, coords);

    // Now coords[0] is the row index, and coords[1] is the column index of the current process

    // Example: Exchange information between neighboring processes in the top and bottom directions
    int top_neighbor, bottom_neighbor, left_neighbor, right_neighbor;
    MPI_Cart_shift(grid_comm, 0, 1, &top_neighbor, &bottom_neighbor);
    MPI_Cart_shift(grid_comm, 1, 1, &left_neighbor, &right_neighbor);

    int upperRightCoords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] + 1) % dims[1]};
    int lowerLeftCoords[2] = {(coords[0] + 1) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};
    int lowerRightCoords[2] = {(coords[0] + 1) % dims[0], (coords[1] + 1) % dims[1]};
    int upperLeftCoords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};

    int upperRight_neighbor, lowerLeft_neighbor, lowerRight_neighbor, upperLeft_neighbor;

    MPI_Cart_rank(grid_comm, upperRightCoords, &upperRight_neighbor);
    MPI_Cart_rank(grid_comm, lowerLeftCoords, &lowerLeft_neighbor);
    MPI_Cart_rank(grid_comm, lowerRightCoords, &lowerRight_neighbor);
    MPI_Cart_rank(grid_comm, upperLeftCoords, &upperLeft_neighbor);

    int neighbors[8] = {top_neighbor, bottom_neighbor, left_neighbor, right_neighbor,
                        upperRight_neighbor, lowerLeft_neighbor, lowerRight_neighbor, upperLeft_neighbor};

    if (rank==4) {
        std::cout << "rank_rows: " << rank_rows << " rank_cols: " << rank_cols << std::endl;
        std::cout << "coords[0]: " << coords[0] << " coords[1]: " << coords[1] << std::endl;
        std::cout << "top_neighbor: " << top_neighbor << " bottom_neighbor: " << bottom_neighbor << std::endl;
        std::cout << "left_neighbor: " << left_neighbor << " right_neighbor: " << right_neighbor << std::endl;
        std::cout << "upperRight_neighbor: " << upperRight_neighbor << " lowerLeft_neighbor: " << lowerLeft_neighbor << std::endl;
        std::cout << "lowerRight_neighbor: " << lowerRight_neighbor << " upperLeft_neighbor: " << upperLeft_neighbor << std::endl;
    
    }

    // Create a communicator using MPI_Dist_graph_create_adjacent
    MPI_Comm new_comm;
    MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,   // 1. old_communicator
                               8,                // 2. indegree
                               neighbors,        // 3. sources (array of input nodes for indegree)
                               MPI_UNWEIGHTED,   // 4. inweights (weights of the input nodes, in this case, MPI_UNWEIGHTED)
                               8,                // 5. outdegree
                               neighbors,        // 6. destinations (array of output nodes for outdegree)
                               MPI_UNWEIGHTED,   // 7. outweights (weights of the output nodes, in this case, MPI_UNWEIGHTED)
                               MPI_INFO_NULL,    // 8. info (hints and optimization information)
                               0,                // 9. reorder (allows the system to reorder the processes, 0 for false, 1 for true)
                               &new_comm);       // 10. new_communicator (output communicator)

    
    int number_of_ghostlayer_elements = 2*(rank_cols + rank_rows) - 4;
    char* send_data = new char[number_of_ghostlayer_elements]; // Example data to be sent

    // can be done in initilization
    int iterate = rank_cols;
    if (rank_rows > rank_cols)
        iterate = rank_rows;
    //

    //write to buffer from world
    for (int i = 0; i<iterate; i++) { 
        if (i<rank_cols) {
            send_data[i] = world[0][i]; // top left, top, top right
            send_data[i+rank_cols] = world[rank_rows-1][i]; // bottom left, bottom, bottom right
        }
        if (i<rank_rows-2) {
            send_data[i+2*rank_cols] = world[i+1][0]; // left
            send_data[i+2*rank_cols+rank_rows-2] = world[i+1][rank_cols-1]; // right
        }
    }

    // order of ghostcells in send_data: [top left, top, top right, bottom left, bottom, bottom right, left, right]

    int number_of_neighbours = 8;
    // Prepare send and receive buffers
    // can be done in initilization
    int send_counts[8] = {1, rank_cols-2, 1, 1, rank_cols-2, 1, rank_rows-2, rank_rows-2};
    int send_displs[8] = {0};
    for (int i=1; i<8;i++) {
        send_displs[i] = send_counts[i-1] + send_displs[i-1];
    }
    //
    

    // Allocate memory for send and receive buffers
    char* recv_buffer = new char[number_of_ghostlayer_elements];
    if (rank==4) {
        std::cout << "send_data:" << " ";
        for (int i = 0; i<number_of_ghostlayer_elements; i++) {
            std::cout << send_data[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "send_counts:" << " ";
        for (int i = 0; i<8; i++) {
            std::cout << send_counts[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "send_displs:" << " ";
        for (int i = 0; i<8; i++) {
            std::cout << send_displs[i] << " ";
        }
        std::cout << std::endl;
    }

    // Perform MPI_Neighbor_alltoallv
    MPI_Neighbor_alltoallv(&send_data, send_counts, send_displs, MPI_CHAR,
                           recv_buffer, send_counts, send_displs, MPI_CHAR, new_comm);


    //write to world from buffer
    for (int i = 0; i<iterate; i++) { 
        if (i<rank_cols) {
            world[0][i] = recv_buffer[i]; // top left, top, top right
            world[rank_rows-1][i] = recv_buffer[i+rank_cols]; // bottom left, bottom, bottom right
        }
        if (i<rank_rows-2) {
            world[i+1][0] = recv_buffer[i+2*rank_cols]; // left
            world[i+1][rank_cols-1] = recv_buffer[i+2*rank_cols+rank_rows]; // right
        }
    }

    // Output received data
    if (rank==4) {
        std::cout << "Rank " << rank << " has neighbors: ";
        for (int i = 0; i < number_of_neighbours; ++i) {
            std::cout << neighbors[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Rank " << rank << " received data: ";
        for (int i = 0; i < number_of_ghostlayer_elements; ++i) {
            std::cout << recv_buffer[i] << " ";
        }
        std::cout << std::endl;
    }

    // Clean up
    delete[] send_data;
    delete[] recv_buffer;

    // Now you can use MPI_Send and MPI_Recv or MPI_Isend and MPI_Irecv to exchange information with the top and bottom neighbors
    MPI_Comm_free(&new_comm);
    MPI_Finalize();
    return 0;
}