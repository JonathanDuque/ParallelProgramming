#include <stdio.h>
#include <string.h>  /* For strlen             */
#include <mpi.h>     /* For MPI functions, etc */
#include <stdlib.h>

const int MAX_STRING = 100;

/*void get_input(int my_rank, int comm_size, double *A, double *x, int *n, int *iters);*/

void print_vector(char *name, double *y, int m, int my_rank);

void matrix_init(double *local_A, int local_n, int n, int my_rank);

void vector_init(double *v, int n, int my_rank);

/* función para generar <size> cantidad de datos aleatorios */
void gen_data(double *array, int size);

int main(int argc, char *argv[]) {
    //char greeting[MAX_STRING];  /* String storing message */
    int comm_sz;               /* Number of processes    */
    int my_rank;               /* My process rank        */
    int n = 0, seed = 0;

    MPI_Init(&argc, &argv);
    //name dimension and seed
    if (argc == 3) {
        n = strtol(argv[1], NULL, 10);
        seed = strtol(argv[2], NULL, 10);
    }
    srand(seed);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int local_n = n / comm_sz;
    double *local_A, *x, *y;

    local_A = malloc(sizeof(double) * local_n * n);
    x = malloc(sizeof(double) * n);
    y = malloc(sizeof(double) * n);

    // Create and initialize matrix and vectors
    matrix_init(local_A, local_n, n, my_rank);
    vector_init(x, n, my_rank);

    //print_vector("A", local_A, local_n*n, my_rank);
    print_vector("x", x, n, my_rank);
    //print_vector("A", local_A, local_n * n);


    //get_input(my_rank, comm_sz, A, x, &n, &iters);

    if (my_rank != 0) {
        //printf("la dimensión %d desde %d:\n", n, my_rank);
    } else {
        //printf("Real dimensión desde 0: %d\n", n);
    }


    MPI_Finalize();

    return 0;
}  /* main */

void matrix_init(double *local_A, int local_n, int n, int my_rank) {
    //int i;
    double *A = NULL;
    //create vector in proc 0
    if (my_rank == 0) {
        A = malloc(sizeof(double) * n * n);
        gen_data(A, n * n);
        //print_vector("A", A, n * n, 0);
        //distribuiting matrix
        MPI_Scatter(A, local_n * n, MPI_DOUBLE, local_A, local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        free(A);
    } else {
        MPI_Scatter(A, local_n * n, MPI_DOUBLE, local_A, local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}


void vector_init( double *lv, int n, int my_rank) {
    //create vector in proc 0
    //double *v = NULL;
    if (my_rank == 0) {
        v = malloc(sizeof(double) * n);
        gen_data(v, n);
        //distribuiting vector


        //MPI_Allgather(v, n, MPI_DOUBLE, lv, n, MPI_DOUBLE, MPI_COMM_WORLD);
        //free(v);
    }

    //MPI_Bcast(v, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

/*void get_input(int my_rank, int comm_size, double *A, double *x, int *n, int *iters) {
    //int dest;

    if (my_rank == 0) {
        int seed;
        printf("Enter n, iterations and seed\n");
        scanf("%d %d %d", n, iters, &seed);

        // la matriz A tendrá una representación unidimensional
        A = malloc(sizeof(double) * *n * *n);
        x = malloc(sizeof(double) * *n);
        //y = malloc(sizeof(double) * n);
        gen_data(A, (*n) * (*n));
        gen_data(x, *n);
    }

    //MPI_Bcast(A, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(x, *n, MPI_DOUBLE, 0, MPI_COMM_WORLD);


}*/

void gen_data(double *array, int size) {
    int i;
    for (i = 0; i < size; i++)
        array[i] = (double) rand() / (double) RAND_MAX;
}

void print_vector(char *name, double *y, int m, int my_rank) {
    int i;
    printf("\nVector %s from %d\n", name, my_rank);
    for (i = 0; i < m; i++)
        printf("%f ", y[i]);
    printf("\n");
}