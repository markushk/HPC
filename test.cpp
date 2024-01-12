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



void print_field(matrix2D<char> &world, int rows, int cols) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << " " << world(i,j);
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

void initial_configuration(matrix2D<char> &world, int rows, int cols, int seed, double probability) {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            world(i,j) = get_random_value(i, j, cols, seed, probability);
        }
    }
}

int c(int x, int y, int cols) {
    return x * cols + y;
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



    int rank_rows = rows/px + 2;
    int rank_cols = cols / py + 2;

    matrix2D<char> world(rank_rows, rank_cols, '0'+rank);
    matrix2D<char> world_copy(rank_rows, rank_cols, '.');
    // initial_configuration(world, rank_rows, rank_cols, seed, probability/100);

    if (rank == 0) {
        std::cout << "initial configuration: \n";
        print_field(world, rank_rows, rank_cols);
    }

    // Create a Cartesian communicator for the 2D grid
    MPI_Comm grid_comm;
    int dims[2] = {px, py};
    int periods[2] = {1, 1}; // Periodic grid in both dimensions
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &grid_comm);

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

    int neighbors[8] = {upperLeft_neighbor, lowerLeft_neighbor, upperRight_neighbor, lowerRight_neighbor,
                        top_neighbor, bottom_neighbor, left_neighbor, right_neighbor};

    if (rank == 0) {
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

    


    int t = 8;
    int sendcount[t], recvcount[t];
    MPI_Aint senddisp[t], recvdisp[t];
    MPI_Datatype sendtype[t], recvtype[t];

    for (int i=0; i<t; i++) {
        sendcount[i] = 1;
        recvcount[i] = 1;
    }

    MPI_Datatype ROW, COL, COR;
  
    MPI_Type_contiguous(rank_cols-2,MPI_CHAR,&ROW);
    MPI_Type_commit(&ROW);
    MPI_Type_vector(rank_rows-2,1,rank_cols,MPI_CHAR,&COL);
    MPI_Type_commit(&COL);
    MPI_Type_contiguous(1,MPI_CHAR,&COR);
    MPI_Type_commit(&COR);

    if (rank == 0) {
        std::cout << "ROW: " << rank_cols-2;
        std::cout << "\nCOL: " << rank_rows-2 << std::endl;
    }


    // order: top left, top right, bottom left, bottom right, left, right, top, bottom,  
    
    senddisp[0] = c(1, 1, rank_cols);                       sendtype[0] = COR; // top left
    senddisp[1] = c(1, rank_cols-2, rank_cols);             sendtype[1] = COR; // top right
    senddisp[2] = c(rank_rows-2, 1, rank_cols);             sendtype[2] = COR; // bottom left
    senddisp[3] = c(rank_rows-2, rank_cols-2, rank_cols);   sendtype[3] = COR; // bottom left
    senddisp[4] = senddisp[0];                              sendtype[4] = COL; // left
    senddisp[5] = senddisp[1];                              sendtype[5] = COL; // right
    senddisp[6] = senddisp[0];                              sendtype[6] = ROW; // top
    senddisp[7] = senddisp[2];                              sendtype[7] = ROW; // bottom

    // ghost layers:
    // order: bottom right, bottom left, top right, top left, right, left, bottom, top,
    recvdisp[0] = c(rank_rows-1, rank_cols-1, rank_cols);   recvtype[0] = COR; // bottom right
    recvdisp[1] = c(rank_rows-1, 0, rank_cols);             recvtype[1] = COR; // bottom left
    recvdisp[2] = c(0, rank_cols-1, rank_cols);             recvtype[2] = COR; // top right
    recvdisp[3] = c(0, 0, rank_cols);                       recvtype[3] = COR; // top left
    recvdisp[4] = c(1, rank_cols-1, rank_cols);             recvtype[4] = COL; // right
    recvdisp[5] = c(1, 0, rank_cols);                       recvtype[5] = COL; // left
    recvdisp[6] = c(rank_rows-1, 1, rank_cols);             recvtype[6] = ROW; // bottom
    recvdisp[7] = c(0, 1, rank_cols);                       recvtype[7] = ROW; // top
        
    // byte offsets
    for (int i=0; i<t; i++) {
        senddisp[i] *= sizeof(char);
        recvdisp[i] *= sizeof(char);
    }

    if (rank == 0) {
        std::cout << "send_data:" << "\n";
        print_field(world, rank_rows, rank_cols);

        std::cout << "sendcount:" << " ";
        for (int i = 0; i<8; i++) {
            std::cout << sendcount[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "senddisp:" << " ";
        for (int i = 0; i<8; i++) {
            std::cout << senddisp[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "recvcount:" << " ";
        for (int i = 0; i<8; i++) {
            std::cout << recvcount[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "recvdisp:" << " ";
        for (int i = 0; i<8; i++) {
            std::cout << recvdisp[i] << " ";
        }
        std::cout << std::endl;
    }


    MPI_Neighbor_alltoallw(&world(0,0),sendcount,senddisp,sendtype,
			   &world(0,0),recvcount,recvdisp,recvtype,new_comm);



    // Output received data
    if (rank == 0) {
        std::cout << "Rank " << rank << " has neighbors: ";
        for (int i = 0; i < t; ++i) {
            std::cout << neighbors[i] << " ";
        }
        std::cout << std::endl;

        std::cout << "Rank " << rank << " received data: \n";
        print_field(world, rank_rows, rank_cols);
        std::cout << std::endl;
    }

    // Clean up

    // Now you can use MPI_Send and MPI_Recv or MPI_Isend and MPI_Irecv to exchange information with the top and bottom neighbors
    MPI_Comm_free(&new_comm);
    MPI_Finalize();
    return 0;
}