//
// Created by Jonathan Duque Gallego on 10/19/19.
//

#include <stdio.h>
#include <stdlib.h>


#include <limits.h>

#include <string.h>
#include <omp.h>
#include <assert.h>

void gen_data(int *array, int size);
void make_sum(int start, int end, int max_limit, int *array, double *sum);

int main(int argc, char *argv[]) {
    int max = 800000000;
    int total_threads = 4;
    int *v;
    int i;
    double sum = 0;
    double time;
    assert(max%total_threads == 0); //299695599691.00 1.8

    v = malloc(sizeof(int) * max);
    srand(1);
    gen_data(v, max);

    time = omp_get_wtime();

#pragma omp parallel num_threads(total_threads) reduction(+:sum)
    make_sum(0, max, max, v, &sum);

    time = omp_get_wtime() - time;
    printf("\nExecution time : %.2f seconds\n", time);

    /*
    for (i = 0; i < max; i++)
        printf("%d ", v[i]);*/

    printf("\nResult: %.2f \n", sum);

    return 0;
}

void gen_data(int *array, int size) {
    int i;
    //#pragma omp parallel for
    for (i = 0; i < size; i++) {
        array[i] = (int) rand() % 500;
        //printf("%d ", array[i]);
    }
}

void make_sum(int start, int end, int max_limit, int *array, double *sum){
    int local_start, local_end;
    int i=0;
    double my_sum=0;
    int my_rank = omp_get_thread_num();
    int total_thread = omp_get_num_threads();
    int partitions = max_limit/total_thread;
    local_start = my_rank*partitions;
    local_end = local_start + partitions;

    for (i = local_start; i < local_end; i++) {
        my_sum += array[i];
    }

    printf("%d My Result: %.2f \n", my_rank, my_sum);

//#pragma omp critical
    *sum += my_sum;
}

