#!/bin/bash
gcc -Wall -fopenmp -o my_it_mat_vect_mult my_it_mat_vect_mult.c
./my_it_mat_vect_mult 5000 100 1 1
./my_it_mat_vect_mult 5000 100 1 2
./my_it_mat_vect_mult 5000 100 1 4
./my_it_mat_vect_mult 5000 100 1 8
./my_it_mat_vect_mult 5000 100 1 16

./my_it_mat_vect_mult 10000 100 1 1
./my_it_mat_vect_mult 10000 100 1 2
./my_it_mat_vect_mult 10000 100 1 4
./my_it_mat_vect_mult 10000 100 1 8

./my_it_mat_vect_mult 15000 100 1 1
./my_it_mat_vect_mult 15000 100 1 2
./my_it_mat_vect_mult 15000 100 1 4
./my_it_mat_vect_mult 15000 100 1 8

./my_it_mat_vect_mult 20000 100 1 1
./my_it_mat_vect_mult 20000 100 1 2
./my_it_mat_vect_mult 20000 100 1 4
./my_it_mat_vect_mult 20000 100 1 8

./my_it_mat_vect_mult 25000 100 1 1
./my_it_mat_vect_mult 25000 100 1 2
./my_it_mat_vect_mult 25000 100 1 4
./my_it_mat_vect_mult 25000 100 1 8
