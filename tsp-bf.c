/**
 *   \file tsp-bf.c
 *   \brief Very simple brute-force algorithm implementation to solve TSP
 *  
 *   This program implements the brute-force algorithm to solve the 
 *   Traveling Salesman Problem. This implementation recursively computes
 *   all possible permutations (Hamiltonian cycles or tour) for a given 
 *   vector. The cost associated to a given permutation is computed and 
 *   the minimum cost is updated to a global variable.
 *   
 *   compile: gcc -Wall -fopenmp -o tsp-bf tsp-bf.c
 *   
 *   \author: Danny Munera (2018)
 *            Parallel Programming course Universidad de Antioquia
 * 
 *   Disclaimer: this implementation was partialy based on this code:
 *  http://wikistack.com/traveling-salesman-problem-brute-force-and-dynamic-programming/ 
 *
 *  TODO: print best tour vector 
 *
 */

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include <math.h>

//#define DEBUG

// global variable containing the minimum cost found
int min_cost = INT_MAX;
// global variable pointing to best tour found
int *best_tour;

void swap(int *x, int *y);

int compute_cost(int *tour, int *dist, int n);

void brute_force_sequential(int *tour, int *dist, int start, int end);

void brute_force_recursive(int *tour, int *dist, int start, int end, int *my_best_tour, int *my_min_cost);

void brute_force_parallel(int *tour, int *dist, int start, int end, int *tasks_assign);

int *get_my_tasks(int my_rank, int *tasks_assign, int maximum_tasks, int size);

int main(int argc, char *argv[]) {
    int size = 4;
    long seed;
    int i, j;
    double time;

    if (argc < 2 && argc > 3) {
        printf("Usage: %s <size> <seed>", argv[0]);
    }

    size = strtol(argv[1], NULL, 10);
    seed = argc == 3 ? strtol(argv[2], NULL, 10) : 1L;

    int iterations = strtol(argv[3], NULL, 10);
    int total_threads = strtol(argv[4], NULL, 10);

    //random number seed
    srand(seed);

    //Distance matrix initialization
    int *dist = (int *) malloc(sizeof(int) * size * size);
    for (i = 0; i < size; i++) {
        for (j = 0; j < size; j++) {
            if (i == j) {
                // accessing matrix throught vector notation
                // dit[i][j] => dist[i*size+j]
                dist[i * size + j] = 0;
            } else {
                dist[i * size + j] = rand() % 10;
            }
        }
    }

#ifdef DEBUG
    printf("Dist matrix:\n");
    for (i = 0; i < size; i++){
      for (j = 0; j < size; j++){
        printf("%d ",dist[i*size+j]);
      }
      printf("\n");
    }
    printf("\n");
#endif // DEBUG

    //best_tour vector initialization
    best_tour = (int *) malloc(sizeof(int) * size);

    //tour vector
    int *tour = (int *) malloc(sizeof(int) * size);
    //tour initialization

    //this array saved the tour changes assigned for each thread
    int *tasks_assigned_for_thread = (int *) malloc(sizeof(int) * size);
    for (i = 0; i < size; i++) {
        tour[i] = i;
        //the module operation gives a number between 0 and total_threads-1
        //this number will be the thread responsible two execute the change that represent the index of
        // tasks_assigned_for_thread array
        tasks_assigned_for_thread[i] = i % total_threads;
    }

    time = omp_get_wtime();
    for (int k = 0; k < iterations ; ++k) {
#pragma omp parallel num_threads(total_threads)
        brute_force_parallel(tour, dist, 0, size - 1, tasks_assigned_for_thread);
    }
    time = omp_get_wtime() - time;

    printf("\nMinimum cost %d\ntour: ", min_cost);
    for (i = 0; i < size; i++)
        printf("%d ", best_tour[i]);
    printf("\nExecution time: %.2f seconds\n", time/iterations);

    return 0;
}

/**
 *  \brief swap
 *  Swap values form two given positions on a vector 
 *  \param x    
 *  \param y 
 */
void swap(int *x, int *y) {
    int temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

/**
 *  \brief compute_cost
 *  Compute cost of a given TSP tour vector
 *  \param tour: pointer to the tour vector
 *  \param dist: pointer to the distance matrix
 *  \param size: size of the tour vector 
 *  \return cost: cost of the tour
 */
int compute_cost(int *tour, int *dist, int size) {
    int i, cost = 0, index;
    for (i = 0; i < size; i++) {
        // distance between cities i and i+1
        // dist[ tour[i%size] ] [ tour[(i+1)%size] ]
        index = tour[i % size] * size + tour[(i + 1) % size];
        cost += dist[index];
    }
#ifdef DEBUG
    printf("Tour ");
    for(i = 0; i < size;i++){
      printf("%d ", tour[i]);
    }
    printf("cost = %d\n",cost);
#endif // DEBUG
    return cost;
}

/**
 *  \brief brute_force
 *  Brute force algorithm to solve the TSP problem
 *  This function recursively computes all possible permutations for
 *  a given TSP instance and updates the minimum cost and best vector to 
 *  global variables (min_cost, best_tour)
 *  \param tour pointer to the tour vector}
 *  \param dist pointer to the distance matrix
 *  \param start tour vector initial position 
 *  \param end tour vector last position
 *
 */

void brute_force_parallel(int *tour, int *dist, int start, int end, int *tasks_assign) {
    //size  = end + 1
    int my_rank = omp_get_thread_num();
    int total_thread = omp_get_num_threads();

    int maximum_task_for_thread = (int) floor((end + 1) / total_thread);

    int *my_tasks;
    //this occur when the load for this thread is unbalanced
    if ((end + 1) % total_thread != 0)
        maximum_task_for_thread += 1;
    //printf(" maximum task %d ", maximum_task_for_thread);

    //getting the array tasks corresponding to my rank
    my_tasks = get_my_tasks(my_rank, tasks_assign, maximum_task_for_thread, end + 1);

    /*
#pragma omp critical
    printf("Thread %d: my tasks.\n", my_rank);
    for (int i = 0; i < maximum_task; i++)
        printf("%d ", my_tasks[i]);
    printf("\n ");
     */

    //every thread has its best_tour , tour and min_cost
    int my_min_cost = INT_MAX;
    int *my_best_tour = (int *) malloc(sizeof(int) * end + 1);
    int *my_tour = (int *) malloc(sizeof(int) * end + 1);
    memcpy(my_tour, tour, sizeof(int) * (end + 1));

    for (int i = 0; i < maximum_task_for_thread; i++) {
        //printf("\nfrom thread %D change %d for %d\n", my_rank, my_tour[start], my_tour[i]);
        if (my_tasks[i] == -1)
            continue;
        swap(&my_tour[start], &my_tour[my_tasks[i]]);
        brute_force_recursive(my_tour, dist, start + 1, end, my_best_tour, &my_min_cost);
        swap(&my_tour[start], &my_tour[my_tasks[i]]); //return change
    }

    int cost = compute_cost(my_best_tour, dist, end + 1);

    //printf("\nThread %d finished\n", my_rank);
    // Compute cost of each permutation
#pragma omp critical
    if (min_cost > cost) {
        // Best solution found - copy cost and tour
        min_cost = cost;
        //update global variable best_tour
        memcpy(best_tour, my_best_tour, sizeof(int) * (end + 1));
    }
}

void brute_force_sequential(int *tour, int *dist, int start, int end) {
    int i, cost;
    if (start == end) {
        // Compute cost of each permutation
        cost = compute_cost(tour, dist, end + 1);
        if (min_cost > cost) {
            // Best solution found - copy cost and tour
            min_cost = cost;
            memcpy(best_tour, tour, sizeof(int) * (end + 1));
        }
    } else {
        for (i = start; i <= end; i++) {
            swap(&tour[start], &tour[i]);
            brute_force_sequential(tour, dist, start + 1, end);
            swap(&tour[start], &tour[i]); //return change
        }
    }
}

void brute_force_recursive(int *tour, int *dist, int start, int end, int *my_best_tour, int *my_min_cost) {
    int i, cost;
    if (start == end) {
        // Compute cost of each permution
        cost = compute_cost(tour, dist, end + 1);

        if (*my_min_cost > cost) {
            // Best solution found - copy cost and tour
            *my_min_cost = cost;
            memcpy(my_best_tour, tour, sizeof(int) * (end + 1));
        }
    } else {
        for (i = start; i <= end; i++) {
            swap(&tour[start], &tour[i]);
            brute_force_recursive(tour, dist, start + 1, end, my_best_tour, my_min_cost);
            swap(&tour[start], &tour[i]); //return change
        }
    }
}

int *get_my_tasks(int my_rank, int *tasks_assign, int maximum_tasks, int size) {
    int *my_tasks = (int *) malloc(sizeof(int) * maximum_tasks);

    for (int j = 0; j < maximum_tasks; ++j) {
        my_tasks[j] = -1;
    }

    int p = 0;
    for (int i = 0; i < size; i++) {
        if (tasks_assign[i] == my_rank) {
            my_tasks[p] = i;
            p++;
        }
    }

    return my_tasks;
}
