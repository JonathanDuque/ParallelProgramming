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
#include <limits.h>
#include <string.h>
#include <omp.h>
#include <assert.h>

/* function to generate <size> amount of random data */
void gen_data(double *array, int size);

void mat_vect_mult(double *A, double *x, double *y, int n, int it);

void mat_vect_mult_parallel(double *A, double *x, double *y, int n, int it);

void print_vector(char *name, double *y, int m);

int main(int argc, char *argv[]) {
    double time;
    double *A = NULL;
    double *x = NULL;
    double *y = NULL;
    int n, iters;
    long seed;
    int total_threads;

    // get dimensions
    n = strtol(argv[1], NULL, 10);
    // get iterations
    iters = strtol(argv[2], NULL, 10);
    // get seed
    seed = strtol(argv[3], NULL, 10);
    srand(seed);

    total_threads = strtol(argv[4], NULL, 10);
    printf("Data parameters: size: %d  iters: %d   seed: %ld   threads: %d", n, iters, seed, total_threads);

    assert(n % total_threads == 0);//check dimension is multiple of number of threads

    //matrix is  unidimensional
    A = malloc(sizeof(double) * n * n);
    x = malloc(sizeof(double) * n);
    y = malloc(sizeof(double) * n);

    //generar valores para las matrices
    time = omp_get_wtime();
    gen_data(A, n * n);
    gen_data(x, n);
    time = omp_get_wtime() - time;
    printf("\nInitialization time  : %.2f seconds\n", time);

    //time = omp_get_wtime();
    //mat_vect_mult(A, x, y, n, iters);
    //time = omp_get_wtime() - time;
    //printf("Execution time sequential : %.2f seconds\n", time);
    //print_vector("y Sequential", y, n);

    double total_time = 0;
    int experiments = 1;

    for (int j = 0; j < experiments; j++) {
        time = omp_get_wtime();
        for (int h = 0; h < iters; h++) {
#pragma omp parallel num_threads(total_threads)
            mat_vect_mult_parallel(A, x, y, n, iters);
//#pragma omp barrier //all thread should be finished
#pragma omp parallel for // no make difference, the processor parallelize automatically!!
//unrolling loop
            for (int i = 0; i < n; i += 4) {
                x[i] = y[i];
                x[i + 1] = y[i + 1];
                x[i + 2] = y[i + 2];
                x[i + 3] = y[i + 3];
            }
        }
        time = omp_get_wtime() - time;
        total_time+=time;
    }
    //printf("Execution time parallel: %.2f seconds\n", time);
    printf("Execution time parallel average: %.2f seconds\n", total_time/experiments);

    //print_vector("y", y, n);
    free(x);
    free(y);

    return 0;
}

void gen_data(double *array, int size) {
    int i;
    //#pragma omp parallel for produce a different random always
    for (i = 0; i < size; i++)
        array[i] = (double) rand() / (double) RAND_MAX;
}

void mat_vect_mult(double *A, double *x, double *y, int n, int it) {
    int h, i, j;
    for (h = 0; h < it; h++) {
        for (i = 0; i < n; i++) {
            y[i] = 0.0;
            //multiplication a A row for vector X and save in y vector position
            for (j = 0; j < n; j++)
                y[i] += A[i * n + j] * x[j];
        }
        // x <= y
        for (i = 0; i < n; i++)
            x[i] = y[i];
    }
}

void mat_vect_mult_parallel(double *A, double *x, double *y, int n, int it) {
    int i, j;
    int local_start, local_end;
    double *my_local_y = NULL;

    int my_rank = omp_get_thread_num();
    int total_thread = omp_get_num_threads();

    int local_n = n / total_thread;

    my_local_y = malloc(sizeof(double) * local_n);
    local_start = my_rank * local_n;
    local_end = local_start + local_n;

    for (i = local_start; i < local_end; i++) {
        my_local_y[i - local_start] = 0.0;
        //multiplication a A row for vector X and save in y vector position
        for (j = 0; j < n; j++)
            my_local_y[i - local_start] += A[i * n + j] * x[j];
    }

#pragma omp critical
    for (i = local_start; i < local_end; i++)
        y[i] = my_local_y[i - local_start];
    // x <= y

    //print_vector(" my y", my_local_y, local_n);

    free(my_local_y);
}

void print_vector(char *name, double *y, int m) {
    int i;
    printf("\nVector %s\n", name);
    for (i = 0; i < m; i++)
        printf("%f ", y[i]);
    printf("\n");
}
