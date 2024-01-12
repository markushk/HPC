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

void stencil_orig(MPI_Comm comm, int rows, int cols, int px, int py)
{
  int d, p, n;
  MPI_Comm cart, cartcomm;
  int rank;
  int i, j;
  int rank_rows = rows/py + 2;
  int rank_cols = cols / px + 2;
  


  d = 2; // number of dimensions

  int order[2] = {px, py};
  int periodic[2] = {1,1};
  int reorder;

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
    stencil_orig(stencomm, rows, cols, px, py);

    MPI_Finalize();
    return 0;
}