#include <stdio.h>
#include <string.h>  /* For strlen             */
#include <mpi.h>     /* For MPI functions, etc */
#include <stdlib.h>

#define DEBUG


/*void get_input(int my_rank, int comm_size, double *A, double *x, int *n, int *iters);*/

void print_vector(char *name, double *y, int n, int my_rank);

void print_vector2(double *local_y, int local_n, int n, char *title, int my_rank);

void matrix_init(double *local_A, int local_n, int n, int my_rank);

void vector_init(double *v, int n, int my_rank);

/* función para generar <size> cantidad de datos aleatorios */
void gen_data(double *array, int size);

void my_bcast(void *data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator);

int main(int argc, char *argv[]) {
    //char greeting[MAX_STRING];  /* String storing message */
    int comm_sz;               /* Number of processes    */
    int my_rank;               /* My process rank        */
    int n = 1, seed = 1;
    int i, j;
    double local_start, local_finish;

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
    double *local_A, *x, *local_y;

    local_A = malloc(sizeof(double) * local_n * n);
    x = malloc(sizeof(double) * n);
    local_y = malloc(sizeof(double) * local_n);
    //printf("Greetings from process %d of %d!\n", my_rank, comm_sz);

    // Create and initialize matrix and vectors
    vector_init(x, n, my_rank);
    matrix_init(local_A, local_n, n, my_rank);
    MPI_Barrier(MPI_COMM_WORLD);
#ifdef DEBUG
    print_vector("A", local_A, local_n*n, my_rank);
    print_vector("x", x, n, my_rank);
#endif // DEBUG

    local_start = MPI_Wtime();
    /*Inicio del código al cual queremos tomarle el tiempo */
    for (i = 0; i < local_n; i++) {
        local_y[i] = 0.0;
        for (j = 0; j < n; j++)
            local_y[i] += local_A[i * n + j] * x[j];
    }
    MPI_Barrier(MPI_COMM_WORLD);//no continue until all process finish
    /*Fin del código al cual queremos tomarle el tiempo*/
    local_finish = MPI_Wtime();

    if (my_rank == 0) {
        // Solo el proceso 0 imprime el tiempo transcurrido
        printf("Tiempo de ejecución = %5.2f segundos \n", (local_finish - local_start));
    }

    print_vector2(local_y, local_n, n, "Y ", my_rank);

    /*double *y;
    if (my_rank == 0) {
        y =  malloc(sizeof(double) * n);
        MPI_Gather(local_y, local_n, MPI_DOUBLE, y, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        printf("%s", title);
        for (i = 0; i < n; i++) {
            printf("%2d ", vector[i]);
        }
        printf("\n");
        free(vector);
    } else {
        MPI_Gather(lv, ln, MPI_INT, vector, ln, MPI_INT, 0, comm);
    }*/




    //get_input(my_rank, comm_sz, A, x, &n, &iters);

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

void vector_init(double *v, int n, int my_rank) {
    //create vector in proc 0
    //double *v = NULL;
    if (my_rank == 0) {
        v = malloc(sizeof(double) * n);
        gen_data(v, n);
        //distribuiting vector

        //MPI_Allgather(v, n, MPI_DOUBLE, lv, n, MPI_DOUBLE, MPI_COMM_WORLD);
        //free(v);
        //my_bcast(v,n,MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Bcast(v, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    }
    MPI_Bcast(v, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //my_bcast(v,n,MPI_DOUBLE, 0, MPI_COMM_WORLD);

}

void my_bcast(void *data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator) {
    int world_rank;
    MPI_Comm_rank(communicator, &world_rank);
    int world_size;
    MPI_Comm_size(communicator, &world_size);

    if (world_rank == root) {
        // If we are the root process, send our data to everyone
        int i;
        for (i = 0; i < world_size; i++) {
            if (i != world_rank) {
                MPI_Send(data, count, datatype, i, 0, communicator);
            }
        }
    } else {
        // If we are a receiver process, receive the data from the root
        MPI_Recv(data, count, datatype, root, 0, communicator,
                 MPI_STATUS_IGNORE);
    }
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

void print_vector(char *name, double *y, int n, int my_rank) {
    int i;
    printf("\nVector %s from %d\n", name, my_rank);
    for (i = 0; i < n; i++)
        printf("%f ", y[i]);
    printf("\n");
}

void print_vector2(double *local_y, int local_n, int n, char *title, int my_rank) {
    int i;
    double *vector = NULL;

    if (my_rank == 0) {
        vector = (double *) malloc(sizeof(double) * n);
        MPI_Gather(local_y, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //print_vector("Y ", vector, n, my_rank);
        free(vector);
    } else {
        MPI_Gather(local_y, local_n, MPI_DOUBLE, vector, local_n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}