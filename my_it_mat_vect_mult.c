/**
 *   \file my_it_mat_vect_mult.c
 *   \brief Multiplica iterativamente un matriz nxn 
 *          por un vector de n posiciones
 *
 *   \author Danny Múnera
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

const int MAX_STRING = 100;
/* función para generar <size> cantidad de datos aleatorios */
void gen_data(double *array, int size);

/* función para multiplicar iterativamente un matriz
 * <m x n> por un vector de tam <n> */
void mat_vect_mult(double *A, double *x, double *y, int n, int it);

/* función para imprimir un vector llamado <name> de tamaño <m>*/
void print_vector(char *name, double *y, int m);

int main(int argc, char *argv[]) {
    double time;
    double *A = NULL;
    double *x = NULL;
    double *y = NULL;
    int n, iters;
    long seed;
    int total_threads;

    // Obtener las dimensiones
    //printf("Ingrese la dimensión n:\n");
    //scanf("%d", &n);
    n = strtol(argv[1], NULL, 10);
    //printf("Ingrese el número de iteraciones:\n");
    //scanf("%d", &iters);
    iters = strtol(argv[2], NULL, 10);
    //printf("Ingrese semilla para el generador de números aleatorios:\n");
    //scanf("%ld", &seed);
    seed = strtol(argv[3], NULL, 10);
    srand(seed);

    total_threads = strtol(argv[4], NULL, 10);
    printf("Data parameters: size: %d  iters: %d   seed: %ld   threads: %d\n", n, iters, seed, total_threads);

    // la matriz A tendrá una representación unidimensional
    A = malloc(sizeof(double) * n * n);
    x = malloc(sizeof(double) * n);
    y = malloc(sizeof(double) * n);

    //generar valores para las matrices
    gen_data(A, n * n);
    gen_data(x, n);

    char       greeting[MAX_STRING];  /* String storing message */
    int        comm_sz;               /* Number of processes    */
    int        my_rank;

    //time = mpi_get_wtime();
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank != 0) {
        sprintf(greeting, "Greetings2 from process %d of %d!", my_rank, comm_sz);
        MPI_Send(greeting, strlen(greeting)+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
    } else {
        printf("Greetings from process %d of %d!\n", my_rank, comm_sz);
        for (int q = 1; q < comm_sz; q++) {
            MPI_Recv(greeting, MAX_STRING, MPI_CHAR, q,0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            printf("%s\n", greeting);
        }
    }
    mat_vect_mult(A, x, y, n, iters);
    MPI_Finalize();

    printf("Execution time parallel average: %.2f seconds\n", 0.0);


    //print_vector("y", y, n);
    free(A);
    free(x);
    free(y);

    return 0;
}

void gen_data(double *array, int size) {
    int i;
    for (i = 0; i < size; i++)
        array[i] = (double) rand() / (double) RAND_MAX;
}

void mat_vect_mult(double *A, double *x, double *y, int n, int it) {
    int h, i, j;
    for (h = 0; h < it; h++) {
        for (i = 0; i < n; i++) {
            y[i] = 0.0;
            for (j = 0; j < n; j++)
                y[i] += A[i * n + j] * x[j];
        }
        // x <= y
        for (i = 0; i < n; i++)
            x[i] = y[i];
    }
}

void print_vector(char *name, double *y, int m) {
    int i;
    printf("\nVector %s\n", name);
    for (i = 0; i < m; i++)
        printf("%f ", y[i]);
    printf("\n");
}
