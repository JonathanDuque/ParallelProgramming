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

/* función para generar <size> cantidad de datos aleatorios */
void gen_data(double *array, int size);

/* función para multiplicar iterativamente un matriz
 * <m x n> por un vector de tam <n> */
void mat_vect_mult(double *A, double *x, double *y, int n, int it);

void mat_vect_mult_parallel(double *A, double *x, double *y, int n, int it);

/* función para imprimir un vector llamado <name> de tamaño <m>*/
void print_vector(char *name, double *y, int m);

int main(int argc, char *argv[]) {
    double time;
    double *A = NULL;
    double *x = NULL;
    double *y = NULL;
    int n, iters;
    long seed;
    int total_threads = 4;

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

    assert(n % total_threads == 0);

    // la matriz A tendrá una representación unidimensional
    A = malloc(sizeof(double) * n * n);
    x = malloc(sizeof(double) * n);
    y = malloc(sizeof(double) * n);

    //generar valores para las matrices
    time = omp_get_wtime();
    gen_data(A, n * n);
    gen_data(x, n);
    time = omp_get_wtime() - time;
    printf("\nInitialization time  : %.2f seconds\n", time);


    time = omp_get_wtime();
    //mat_vect_mult(A, x, y, n, iters);
    time = omp_get_wtime() - time;
    printf("Execution time sequential : %.2f seconds\n", time);
    //print_vector("y Sequential", y, n);

    time = omp_get_wtime();

    for (int h = 0; h < iters; h++) {
#pragma omp parallel num_threads(total_threads)
        mat_vect_mult_parallel(A, x, y, n, iters);
#pragma omp parallel for // no make difference, the processor parallelized automatically!!
        for (int i = 0; i < n; i++)
            x[i] = y[i];
    }
    time = omp_get_wtime() - time;
    printf("Execution time parallel: %.2f seconds\n", time); //3973.445822 3952.866606 3998.713330


    // 2 iters     31823826.704084 31860559.788320
    //print_vector("y", y, n);  //2.057209 2.516496 3.516359 3.734857 1.915111 2.585758 2.906468 3.581800 2.101284 3.334560 2.235218 2.519298
    free(A); // 2 iters 15.565738 13.449727 17.000071 21.118892 15.146032 14.727463 17.298086 20.556885 13.282327 19.676767 12.972550 16.944703
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

    //for (i = 0; i < n; i++)
    //x[i] = y[i];

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
