
CPP = /usr/local/opt/llvm/bin/clang
CPPFLAGS = -I/usr/local/opt/llvm/include -fopenmp
LDFLAGS = -L/usr/local/opt/llvm/lib

tsp-bf: tsp-bf.c
#test: test.c
#my_it_mat_vect_mult: my_it_mat_vect_mult.c
	$(CPP) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)
