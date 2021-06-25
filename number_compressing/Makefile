#CC = icc
#CC = g++-4.5
#CC = g++
CC = gcc
FLG = -O3
NAME = pagerank_pthreads

pagerank_pthreads: pagerank_pthreads.c time_manager.c function.c function.h time_manager.h pagerank_pthreads.h

	$(CC) pagerank_pthreads.c time_manager.c function.c  -lpthread -lm -fopenmp -g -o $(NAME) $(FLG)

clean:
	rm -f pagerank_pthreads

run_gdb:
	gdb ./pagerank_pthreads
