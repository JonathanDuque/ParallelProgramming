#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define DEBUG

void vector_init(int *lv, int ln, int n, int my_rank, MPI_Comm comm);

void print_vector(int *lv, int ln, int n, char *title, int rank, MPI_Comm comm);

void vector_sum(int *v1, int *v2, int *res, int n);

int main(int argc, char *argv[]) {
    long seed = 1L;
    int n = 16;
    MPI_Comm comm;
    int ln;
    int *la, *lb, *lr;
    int comm_sz, rank;

    MPI_Init(&argc, &argv);
    if (argc == 3) {
        n = strtol(argv[1], NULL, 10);
        seed = strtol(argv[2], NULL, 10);
    }
    srand(seed);

    comm = MPI_COMM_WORLD;
    MPI_Comm_size(comm, &comm_sz);
    MPI_Comm_rank(comm, &rank);

    ln = n / comm_sz;

    la = (int *) malloc(sizeof(int) * ln);
    lb = (int *) malloc(sizeof(int) * ln);
    lr = (int *) malloc(sizeof(int) * ln);

    // Create and initialize vectors
    vector_init(la, ln, n, rank, comm);
    vector_init(lb, ln, n, rank, comm);
#ifdef DEBUG
    print_vector(la, ln, n, "A ", rank, comm);
    print_vector(lb, ln, n, "B ", rank, comm);
#endif // DEBUG 


    vector_sum(la, lb, lr, ln);

#ifdef DEBUG
    print_vector(lr, ln, n, "R ", rank, comm);
#endif // DEBUG

    MPI_Finalize();

    return 0;
}

void print_vector(int *lv, int ln, int n, char *title, int rank, MPI_Comm comm) {
    int i;
    int *vector = NULL;

    if (rank == 0) {
        vector = (int *) malloc(sizeof(int) * n);
        MPI_Gather(lv, ln, MPI_INT, vector, ln, MPI_INT, 0, comm);
        printf("%s", title);
        for (i = 0; i < n; i++) {
            printf("%2d ", vector[i]);
        }
        printf("\n");
        free(vector);
    } else {
        MPI_Gather(lv, ln, MPI_INT, vector, ln, MPI_INT, 0, comm);
    }
}

void vector_init(int *lv, int ln, int n, int my_rank, MPI_Comm comm) {
    int i;
    int *v = NULL;
    //create vector in proc 0
    if (my_rank == 0) {
        v = (int *) malloc(sizeof(int) * n);
        for (i = 0; i < n; i++) {
            v[i] = rand() % 10;
        }
        //distribuiting vector
        MPI_Scatter(v, ln, MPI_INT, lv, ln, MPI_INT, 0, comm);
        free(v);
    } else {
        MPI_Scatter(v, ln, MPI_INT, lv, ln, MPI_INT, 0, comm);
    }
}

void vector_sum(int *v1,  /*in*/
                int *v2,  /*in*/
                int *res, /*out*/
                int n     /*in*/ ) {
    int i;
    for (i = 0; i < n; i++) {
        res[i] = v1[i] + v2[i];
    }
}
