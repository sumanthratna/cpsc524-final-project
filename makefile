# CC = gcc
# FLG = -O4
# NAME = pagerank_ags

# pagerank_ags: pagerank_parallel_gaussseidel.c pagerank_pthreads.h

# 	$(CC) pagerank_parallel_gaussseidel.c -lpthread -lm -o $(NAME)

# clean:
# 	rm -f *.o *.out *.exe
# 	rm -f *.bin

# CC = gcc
# FLG = -O4
# NAME = pagerank_pthreads

# pagerank_pthreads: pagerank_parallel_jacobi_blocking.c pagerank_pthreads.h

# 	$(CC) pagerank_parallel_jacobi_blocking.c -lpthread -lm -o $(NAME)

# clean:
# 	rm -f *.o *.out *.exe
# 	rm -f *.bin

CC = gcc
FLG = -O4
NAME = pagerank_serial

pagerank_serial: pagerank_serial.cpp pagerank_pthreads.h

	$(CC) pagerank_serial.cpp -lpthread -lm -o $(NAME)

clean:
	rm -f *.o *.out *.exe
	rm -f *.bin