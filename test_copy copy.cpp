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

void stencil_orig(MPI_Comm comm)
{
  int d, p, n;
  MPI_Comm cart, cartcomm;
  int rank;
  int i, j;
  
  MPI_Comm_size(comm,&p);

  d = 2; // number of dimensions

  int order[d];
  int periodic[d];
  int reorder;

  for (i=0; i<d; i++) {
    order[i]    = 0;
    periodic[i] = 1;
  }
  MPI_Dims_create(p,d,order);
  reorder = 0; // reorder here?
  MPI_Cart_create(comm,d,order,periodic,reorder,&cart);
  MPI_Comm_rank(cart,&rank); // reorder could have taken place

  n = 4; // arbitrary choice of local matrix size

  char matrix[n+2][n+2];
  char matrix2[n+2][n+2];

  int t = 8;
  int target[] = { 0,1, 0,-1, -1,0,   1,0,
		  -1,1, 1,1,  1,-1, -1,-1};

  int sources[t], destinations[t];
  for (i=0; i<t; i++) {
    int vector[d];
    MPI_Cart_coords(cart,rank,d,vector);
    vector[0] += target[d*i];
    vector[1] += target[d*i+1];
    MPI_Cart_rank(cart,vector,&destinations[i]);
    MPI_Cart_coords(cart,rank,d,vector);
    vector[0] -= target[d*i];
    vector[1] -= target[d*i+1];
    MPI_Cart_rank(cart,vector,&sources[i]);
  }  

  if (rank == 0) {
    std::cout << "destinations:" << "\n";
    for (i=0;i<t;i++) {
      std::cout << destinations[i] << " ";
    }
    std::cout << std::endl;
  }

#ifdef DEBUG
  MPI_Comm_rank(comm,&rank);
  if (rank==2) {
    printf("Prop: ");
    for (i=0; i<t; i++) printf("%d ",destinations[i]);
    printf("\n");
  }
#endif
  
  reorder = 0; // or reorder here?     
  MPI_Dist_graph_create_adjacent(cart,
				 t,destinations,MPI_UNWEIGHTED,
				 t,destinations,MPI_UNWEIGHTED,
				 MPI_INFO_NULL,reorder,&cartcomm);

  MPI_Datatype ROW, COL, COR;
  
  MPI_Type_contiguous(n,MPI_CHAR,&ROW);
  MPI_Type_commit(&ROW);
  MPI_Type_vector(n,1,n+2,MPI_CHAR,&COL);
  MPI_Type_commit(&COL);
  MPI_Type_contiguous(1,MPI_CHAR,&COR);
  MPI_Type_commit(&COR);

  int sendcount[t], recvcount[t];
  MPI_Aint senddisp[t], recvdisp[t];
  MPI_Datatype sendtype[t], recvtype[t];

  // by using datatypes, all counts 1
  for (i=0; i<t; i++) {
    sendcount[i] = 1;
    recvcount[i] = 1;
  }
  // order: right, left, top, bottom, bottom left, top right, bottom right, top left
  senddisp[0] = 1*(n+2)+n;     sendtype[0] = COL;
  recvdisp[0] = 1*(n+2)+n+1;   recvtype[0] = COL;
  senddisp[1] = 1*(n+2)+1;     sendtype[1] = COL;
  recvdisp[1] = 1*(n+2);       recvtype[1] = COL;
  senddisp[2] = 1*(n+2)+1;     sendtype[2] = ROW;
  recvdisp[2] = 1;             recvtype[2] = ROW;
  senddisp[3] = n*(n+2)+1;     sendtype[3] = ROW;
  recvdisp[3] = (n+1)*(n+2)+1; recvtype[3] = ROW;

  senddisp[4] = n*(n+2)+1;       sendtype[4] = COR;
  recvdisp[4] = (n+1)*(n+2);     recvtype[4] = COR;
  senddisp[5] = 1*(n+2)+n;       sendtype[5] = COR;
  recvdisp[5] = n+1;             recvtype[5] = COR;
  senddisp[6] = n*(n+2)+n;       sendtype[6] = COR;
  recvdisp[6] = (n+1)*(n+2)+n+1; recvtype[6] = COR;
  senddisp[7] = 1*(n+2)+1;       sendtype[7] = COR;
  recvdisp[7] = 0;               recvtype[7] = COR;
    
  // byte offsets
  for (i=0; i<t; i++) {
    senddisp[i] *= sizeof(char);
    recvdisp[i] *= sizeof(char);
  }
  

    for (i=0; i<n+2; i++) {
        for (j=0; j<n+2; j++) {
          matrix[i][j] = '0'+rank; // something more sensible
          matrix2[i][j] = '0'+rank;
        }
    }

    if (rank == 0) {
      std::cout << "inital: \n";
      for (int i = 0; i < n+2; ++i) {
          for (int j = 0; j < n+2; ++j) {
              std::cout << " " << matrix2[i][j];
          }
          std::cout << std::endl;
      }
    }

    if (rank == 0) {
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

    // update
    MPI_Neighbor_alltoallw(matrix,sendcount,senddisp,sendtype,
                matrix2,recvcount,recvdisp,recvtype,cartcomm);
        
    if (rank == 0) {
      std::cout << "final: \n";
      for (int i = 0; i < n+2; ++i) {
          for (int j = 0; j < n+2; ++j) {
              std::cout << " " << matrix2[i][j];
          }
          std::cout << std::endl;
      }
    }

  MPI_Comm_free(&cart);
  MPI_Comm_free(&cartcomm);

  MPI_Type_free(&ROW);
  MPI_Type_free(&COL);
  MPI_Type_free(&COR);
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

    matrix2D<char> world(rank_rows, rank_cols, '0'+rank);
    matrix2D<char> world_copy(rank_rows, rank_cols, '.');
    // initial_configuration(world, rank_rows, rank_cols, seed, probability/100);
    MPI_Comm stencomm;
    MPI_Comm_dup(MPI_COMM_WORLD,&stencomm);
    stencil_orig(stencomm);


    // if (rank == 0) {
    //     std::cout << "initial configuration: \n";
    //     print_field(world, rank_rows, rank_cols);
    // }

    // // Create a Cartesian communicator for the 2D grid
    // MPI_Comm grid_comm;
    // int dims[2] = {px, py};
    // int periods[2] = {1, 1}; // Periodic grid in both dimensions
    // MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &grid_comm);

    // // Determine the coordinates of the current process in the 2D grid
    // int coords[2];
    // MPI_Cart_coords(grid_comm, rank, 2, coords);

    // // Now coords[0] is the row index, and coords[1] is the column index of the current process

    // // Example: Exchange information between neighboring processes in the top and bottom directions
    // int top_neighbor, bottom_neighbor, left_neighbor, right_neighbor;
    // MPI_Cart_shift(grid_comm, 0, 1, &top_neighbor, &bottom_neighbor);
    // MPI_Cart_shift(grid_comm, 1, 1, &left_neighbor, &right_neighbor);

    // int upperRightCoords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] + 1) % dims[1]};
    // int lowerLeftCoords[2] = {(coords[0] + 1) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};
    // int lowerRightCoords[2] = {(coords[0] + 1) % dims[0], (coords[1] + 1) % dims[1]};
    // int upperLeftCoords[2] = {(coords[0] - 1 + dims[0]) % dims[0], (coords[1] - 1 + dims[1]) % dims[1]};

    // int upperRight_neighbor, lowerLeft_neighbor, lowerRight_neighbor, upperLeft_neighbor;

    // MPI_Cart_rank(grid_comm, upperRightCoords, &upperRight_neighbor);
    // MPI_Cart_rank(grid_comm, lowerLeftCoords, &lowerLeft_neighbor);
    // MPI_Cart_rank(grid_comm, lowerRightCoords, &lowerRight_neighbor);
    // MPI_Cart_rank(grid_comm, upperLeftCoords, &upperLeft_neighbor);

    // int neighbors[8] = {top_neighbor, bottom_neighbor, left_neighbor, right_neighbor,
    //                     upperRight_neighbor, lowerLeft_neighbor, lowerRight_neighbor, upperLeft_neighbor};

    // if (rank == 0) {
    //     std::cout << "rank_rows: " << rank_rows << " rank_cols: " << rank_cols << std::endl;
    //     std::cout << "coords[0]: " << coords[0] << " coords[1]: " << coords[1] << std::endl;
    //     std::cout << "top_neighbor: " << top_neighbor << " bottom_neighbor: " << bottom_neighbor << std::endl;
    //     std::cout << "left_neighbor: " << left_neighbor << " right_neighbor: " << right_neighbor << std::endl;
    //     std::cout << "upperRight_neighbor: " << upperRight_neighbor << " lowerLeft_neighbor: " << lowerLeft_neighbor << std::endl;
    //     std::cout << "lowerRight_neighbor: " << lowerRight_neighbor << " upperLeft_neighbor: " << upperLeft_neighbor << std::endl;
    
    // }

    // // Create a communicator using MPI_Dist_graph_create_adjacent
    // MPI_Comm new_comm;
    // MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD,   // 1. old_communicator
    //                            8,                // 2. indegree
    //                            neighbors,        // 3. sources (array of input nodes for indegree)
    //                            MPI_UNWEIGHTED,   // 4. inweights (weights of the input nodes, in this case, MPI_UNWEIGHTED)
    //                            8,                // 5. outdegree
    //                            neighbors,        // 6. destinations (array of output nodes for outdegree)
    //                            MPI_UNWEIGHTED,   // 7. outweights (weights of the output nodes, in this case, MPI_UNWEIGHTED)
    //                            MPI_INFO_NULL,    // 8. info (hints and optimization information)
    //                            0,                // 9. reorder (allows the system to reorder the processes, 0 for false, 1 for true)
    //                            &new_comm);       // 10. new_communicator (output communicator)

    


    // int t = 8;
    // int sendcount[t], recvcount[t];
    // MPI_Aint senddisp[t], recvdisp[t];
    // MPI_Datatype sendtype[t], recvtype[t];

    // for (int i=0; i<t; i++) {
    //     sendcount[i] = 1;
    //     recvcount[i] = 1;
    // }

    // MPI_Datatype ROW, COL, COR;
  
    // MPI_Type_contiguous(rank_cols-2,MPI_CHAR,&ROW);
    // MPI_Type_commit(&ROW);
    // MPI_Type_vector(rank_rows-2,1,rank_cols,MPI_CHAR,&COL);
    // MPI_Type_commit(&COL);
    // MPI_Type_contiguous(1,MPI_CHAR,&COR);
    // MPI_Type_commit(&COR);

    // //          0               1           2               3           4           5               6           7
    // // send = [top left,        top,        top right,      left,       right,      bottom left,    bottom,     bottom right]
    // // recv = [bottom right,    bottom,     bottom left,    right,      left,       top right,      top,        top left]
    // //          7               6           5               4           3           2               1           0
    // // -->  recvdisp[7] = senddisp[0]
    // //      recvdisp[6] = senddisp[1]
    // //      ...


    // senddisp[0] = 0;                        sendtype[0] = COR; //top left corner
    // senddisp[7] = rank_cols*rank_rows-1;    sendtype[7] = COR; //bottom right corner
    // recvdisp[0] = senddisp[7];              recvtype[0] = COR; //bottom right corner -> top left corner
    // recvdisp[7] = senddisp[0];              recvtype[7] = COR; //top left corner -> bottom right corner

    // senddisp[2] = rank_cols-1;              sendtype[2] = COR; //top right corner
    // senddisp[5] = rank_cols*(rank_rows-1);  sendtype[5] = COR; //bottom left corner
    // recvdisp[2] = senddisp[5];              recvtype[2] = COR; //bottom left corner -> top right corner
    // recvdisp[5] = senddisp[2];              recvtype[5] = COR; //top right corner -> bottom left corner

    // senddisp[1] = senddisp[0] + 1;          sendtype[1] = ROW; //top
    // senddisp[6] = senddisp[5] + 1;          sendtype[7] = ROW; //bottom
    // recvdisp[1] = senddisp[6];              recvtype[1] = ROW; //bottom -> top
    // recvdisp[6] = senddisp[1];              recvtype[6] = ROW; //top -> bottom

    // senddisp[3] = senddisp[2] + 1;          sendtype[3] = COL; //left
    // senddisp[4] = senddisp[3]+senddisp[2];  sendtype[4] = COL; //right
    // recvdisp[3] = senddisp[4];              recvtype[3] = COL; //right -> left
    // recvdisp[4] = senddisp[3];              recvtype[4] = COL; //left -> right
        
    // // byte offsets
    // for (int i=0; i<t; i++) {
    //     senddisp[i] *= sizeof(char);
    //     recvdisp[i] *= sizeof(char);
    // }

    // if (rank == 0) {
    //     std::cout << "send_data:" << "\n";
    //     print_field(world, rank_rows, rank_cols);

    //     std::cout << "sendcount:" << " ";
    //     for (int i = 0; i<8; i++) {
    //         std::cout << sendcount[i] << " ";
    //     }
    //     std::cout << std::endl;

    //     std::cout << "senddisp:" << " ";
    //     for (int i = 0; i<8; i++) {
    //         std::cout << senddisp[i] << " ";
    //     }
    //     std::cout << std::endl;

    //     std::cout << "recvcount:" << " ";
    //     for (int i = 0; i<8; i++) {
    //         std::cout << recvcount[i] << " ";
    //     }
    //     std::cout << std::endl;

    //     std::cout << "recvdisp:" << " ";
    //     for (int i = 0; i<8; i++) {
    //         std::cout << recvdisp[i] << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // char matrix[5][5];
    // for (int i = 0; i < 5; ++i) {
    //     for (int j = 0; j < 5; ++j) {
    //         matrix[i][j] = '.';
    //     }
    // }

    // MPI_Neighbor_alltoallw(matrix,sendcount,senddisp,sendtype,
	// 		   matrix,recvcount,recvdisp,recvtype,new_comm);



    // // Output received data
    // if (rank == 0) {
    //     std::cout << "Rank " << rank << " has neighbors: ";
    //     for (int i = 0; i < t; ++i) {
    //         std::cout << neighbors[i] << " ";
    //     }
    //     std::cout << std::endl;

    //     std::cout << "Rank " << rank << " received data: \n";
    //     for (int i = 0; i < t; ++i) {
    //         print_field(world_copy, rank_rows, rank_cols);
    //     }
    // }

    // // Clean up

    // // Now you can use MPI_Send and MPI_Recv or MPI_Isend and MPI_Irecv to exchange information with the top and bottom neighbors
    // MPI_Comm_free(&new_comm);
    MPI_Finalize();
    return 0;
}