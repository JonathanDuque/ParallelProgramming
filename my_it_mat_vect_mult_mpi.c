#include <stdio.h>
#include <string.h>  /* For strlen             */
#include <mpi.h>     /* For MPI functions, etc */
#include <stdlib.h>
#include <assert.h>

//#define DEBUG

void print_vector(char *name, double *y, int n, int my_rank);

void print_result(double *local_y, int local_n, int n, char *title, int my_rank);

void matrix_init(double *local_A, int local_n, int n, int my_rank);

void vector_init(double *v, int n, int my_rank);

/* function to generate <size> amount of random data */
void gen_data(double *array, int size);

int main(int argc, char *argv[]) {
    int comm_sz;               /* Number of processes    */
    int my_rank;               /* My process rank        */
    int n = 1, seed = 1, iters = 1;
    int i, j, k;
    double local_start, local_finish, local_elapsed, elapsed;

    MPI_Init(&argc, &argv);

    if (argc == 4) {
        n = strtol(argv[1], NULL, 10); //dimension
        iters = strtol(argv[2], NULL, 10); //iterations number
        seed = strtol(argv[3], NULL, 10); //seed for random numbers
    }
    srand(seed);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    int local_n = n / comm_sz; //used for the local matrix dimension,  matrix processed by each processor
    double *local_A, *x, *local_y;

    assert(n % comm_sz == 0);

    //to ensure all processes init at the same time
    MPI_Barrier(MPI_COMM_WORLD);
    //take init variable time
    local_start = MPI_Wtime();
    local_A = malloc(sizeof(double) * local_n * n);
    x = malloc(sizeof(double) * n);
    local_y = malloc(sizeof(double) * local_n);

    // Create and initialize local matrix and vector
    vector_init(x, n, my_rank);
    matrix_init(local_A, local_n, n, my_rank);

    // finish init variable time
    local_finish = MPI_Wtime();
    local_elapsed = local_finish - local_start;  // each process take a local time
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // take the slow process time

    // process 0 print init time
    if (my_rank == 0) {
        printf("Initialization time = %5.2f seconds \n", elapsed );
    }

#ifdef DEBUG
    print_vector("local A", local_A, local_n*n, my_rank);
    print_vector("x", x, n, my_rank);
#endif // DEBUG

    //to ensure all processes init at the same time
    MPI_Barrier(MPI_COMM_WORLD);
    //take the moment when begin the multiplication
    local_start = MPI_Wtime();
    for (k = 0; k < iters; k++) {
        for (i = 0; i < local_n; i++) {
            local_y[i] = 0.0;
            for (j = 0; j < n; j++)
                local_y[i] += local_A[i * n + j] * x[j];
        }

        // update x value in all process, x=y;
        MPI_Barrier(MPI_COMM_WORLD); //wait until all process finish the multiplication
        //join local y result and put its value in x
        MPI_Gather(local_y, local_n, MPI_DOUBLE, x, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(x, n, MPI_DOUBLE, 0, MPI_COMM_WORLD); //send x to everybody for begin a new iteration
    }
    free(local_A);
    free(x);

    /*take the moment when multiplication end*/
    local_finish = MPI_Wtime();
    local_elapsed = local_finish - local_start;  // each process take a local time
    MPI_Reduce(&local_elapsed, &elapsed, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD); // take the slow process time

    // process 0 print execution time
    if (my_rank == 0) {
        printf("Execution time = %5.2f seconds \n", elapsed);
    }

    print_result(local_y, local_n, n, "Y", my_rank);
    free(local_y);

    MPI_Finalize();

    return 0;
}  /* main */

void matrix_init(double *local_A, int local_n, int n, int my_rank) {
    //int i;
    double *A = NULL;
    //create matrix in proc 0
    if (my_rank == 0) {
        A = malloc(sizeof(double) * n * n);
        gen_data(A, n * n);
        //distribuiting matrix to all process
        MPI_Scatter(A, local_n * n, MPI_DOUBLE, local_A, local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        free(A);
    } else {
        MPI_Scatter(A, local_n * n, MPI_DOUBLE, local_A, local_n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}

void vector_init(double *v, int n, int my_rank) {
    if (my_rank == 0) {
        //fiil vector in process 0
        gen_data(v, n);
    }
    //distribuiting vector
    MPI_Bcast(v, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);//no continue until all process finish
}

void gen_data(double *array, int size) {
    int i;
    for (i = 0; i < size; i++)
        array[i] = (double) rand() / (double) RAND_MAX;
}

void print_vector(char *name, double *y, int n, int my_rank) {
    int i;
    printf("\nVector %s from %d\n", name, my_rank);
    for (i = 0; i < n; i++)
        printf("%f ", y[i]);
    printf("\n");
}

void print_result(double *local_y, int local_n, int n, char *title, int my_rank) {
    int i;
    double *vector = NULL;

    if (my_rank == 0) {
        vector = (double *) malloc(sizeof(double) * n);
        MPI_Gather(local_y, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#ifdef DEBUG
        print_vector(title, vector, n, my_rank);
#endif
        free(vector);
    } else {
        MPI_Gather(local_y, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}