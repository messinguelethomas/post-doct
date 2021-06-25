#ifndef PAGERANK_PTHREADS 
#define PAGERANK_PTHREADS

#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>
#include"time_manager.h"
#include"function.h"
#include <x86intrin.h>

#define taskloop_  -1  
#define pthread_  1  
#define static_   2  
#define dynamic_  3
#define guided_   4
#define auto_     5
#define task_     6
#define pthread_2  7   
#define static_2   8  
#define dynamic_2  9
#define guided_2   10
#define auto_2     11
#define task_2     12
#define taskloop_2  13  

/******************** Structs ********************/
/***** Struct for timestamps ***
struct timeval start,end;*/

/***** Struct used for Threads data *****/

typedef struct
{
	int tid;
	int start, end;
	unsigned long long degree_sum;
} Thread; 

/***** Struct used for Nodes data *****/

typedef unsigned char con_size_t[2];
typedef unsigned char from_size_t[2];
typedef unsigned char node_size_3_t[3];

typedef struct
{
	char size;
	char* node_size;
}node_size_x_t;

typedef struct
{
	double p_t0;
	double p_t1;
	double e;
	int *From_id;
	int con_size;
	int from_size;
}Node;

typedef struct
{
	double p_t0;
	double p_t1;
	double e;
	node_size_3_t *From_id;
	con_size_t con_size;
	from_size_t from_size;
}Node_succint_global;

typedef struct
{
	double p_t0;
	double p_t1;
	double e;
	node_size_x_t *From_id;
	con_size_t con_size;
	from_size_t from_size;
}Node_succint_local;

typedef struct thread_load_info
{
	int tid;
	unsigned long long nb_node_sum;
	unsigned long long degree_sum;
	unsigned long long time;
}thread_load_info_t;

#endif
