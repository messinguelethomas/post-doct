/********************************************************************/
/*    Pagerank project 2014 - Parallel version                      */
/*    	*based on Cleve Moler's matlab implementation               */
/*                                                                  */
/*    Implemented by Nikos Katirtzis (nikos912000)                  */
/********************************************************************/
/*    Modified by thomas Messi Nguélé (messinguelethomas@gmail.com) */

/******************** Includes - Defines ****************/
#include "pagerank_pthreads.h"
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

#define maxThreads 64
#define S_LOCAL 2
#define S_GLOBAL 1
#define S_NONE 0


/******************** Defines ****************/
// Number of nodes
int N, num_threads, real_nb_task, DS;
unsigned long long degree_sum_from = 0, degree_sum_to = 0;
int to_degree = 0, chunk = 1, scd;
tab_load_info_t tab_load, tab_load_reinit, tab_load_lmax,tab_load_ppgr;
// Convergence threashold and algorithm's parameter d  
double threshold, d;

//Table of threads
pthread_t *Threads;

// Table with thread's data
Thread *Threads_data;

// Table of node's data
Node *Nodes;
//begin succint
Node_succint_global*  Nodes_sg;
Node_succint_local*  Nodes_sl;
//end succint
pthread_mutex_t lockP = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t locksum = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lockmax = PTHREAD_MUTEX_INITIALIZER;

// Number of iterations
int iterations = 0;

double max_error = 1;
double sum = 0;

/***** Memory allocation - Initializations for Threads *****/
tab_load_info_t init_tab_load_info(int nb_thread)
{
	int i;
	tab_load_info_t tab_load = (tab_load_info_t)malloc(sizeof(thread_load_info_t)*nb_thread);

	for(i=0; i<nb_thread; i++)
	{
		tab_load[i].nb_node_sum = 0;
		tab_load[i].degree_sum  = 0;
		tab_load[i].time        = 0;
	}
	return tab_load;
}

void write_load_info(char*nom_fich_load, int nb_thread, tab_load_info_t tab_load, int chunk, int scd)
{
	int tid,  i;
        unsigned long long degree = 0, time_sum = 0, degree_sum = 0, temps, time_up = 0, nb_node_sum = 0, nb_node_tid; 
	FILE* time_file = fopen(nom_fich_load, "a");


	for(i=0; i<nb_thread; i++)
	{
		nb_node_tid = tab_load[i].nb_node_sum;
		degree      = tab_load[i].degree_sum;
		temps       = tab_load[i].time      ;
		tid         = i;
		//printf("%d\t%d\t%d\t%lld\t%lld\t%lld.%03ld\n",scd, chunk, tid, nb_node_tid/iterations, degree/iterations, temps/1000,temps%1000);
		fprintf(time_file,"%d\t%d\t%d\t%lld\t%lld\t%lld.%03ld\n",scd, chunk, tid, nb_node_tid/iterations, degree/iterations, temps/1000,temps%1000);
		degree_sum  += degree;
		time_sum    += temps;
		nb_node_sum += nb_node_tid;
	}
	fprintf(time_file,"%d\t%d\t-1\t%lld\t%lld\t%lld.%03ld\n\n",scd, chunk, nb_node_sum/iterations, degree_sum/iterations, time_sum/1000,time_sum%1000);
	//printf("%d\t%d\t-1\t%lld\t%lld\t%lld.%03ld\n\n",scd, chunk, nb_node_sum/iterations, degree_sum/iterations, time_sum/1000,time_sum%1000);
	fclose(time_file);
}

Thread* generate_community_task(char* file_name)
{
	int nb_node = N, num_comm,node_deb,node_fin,nb_node_comm, nb_comm, nb_task;
	unsigned long long degree_sum, sum_degree;
	unsigned long long max_degree_sum_task, degree_sum_temp, i, task_id;
	Thread* task_info;

        FILE* graph_file = NULL,*rslt_file = NULL;
	char* resl_split[100], *nom_fich, *src_temp, *dst_temp, file_name_cpy[100];
        graph_file = fopen(file_name, "r");
	
	fscanf(graph_file, "%d\n", &nb_comm);
	nb_task = nb_comm;
	task_info = (Thread*)malloc(sizeof(Thread)*nb_task);
	
	i = 0; task_id =0;
	while(fscanf(graph_file, "%d\t%d\t%d\t%d\t%lld\n", &num_comm,&node_deb,&node_fin,&nb_node_comm,&sum_degree) != EOF) {
		if(task_id < nb_task)
		{
			task_info[task_id].start = node_deb;
			task_info[task_id].end   = node_fin + 1;
			task_info[task_id].tid = task_id;
			task_info[task_id].degree_sum = sum_degree;
//printf("task %d begin at %d end at %d degree_sum = %lld\n",task_id, task_info[task_id].start, task_info[task_id].end, task_info[task_id].degree_sum);
			task_id++;
			real_nb_task = task_id;
		}
	}

	fclose(graph_file);

	return task_info;
}

Thread* generate_task_index(int to_degree, int nb_thread, int chunk, int ds)
{
	int nb_node = N;
	unsigned long long degree_sum;
	unsigned long long max_degree_sum_task, degree_sum_temp, i, task_id;
	Thread* task_info;
	int nb_task = chunk*nb_thread;

	switch(to_degree)
	{
		case 1: degree_sum = degree_sum_to;
			break;
		case 0: degree_sum = degree_sum_from;
			break;
		default:
			break;
	}

	max_degree_sum_task = degree_sum / nb_task;
	task_info = (Thread*)malloc(sizeof(Thread)*nb_task);

		//printf("degree_sum = %lld max_degree_sum_task = %lld\n",degree_sum, max_degree_sum_task);
	i = 0; task_id =0;
	while((i < nb_node)&&(task_id < nb_task))
	{
		task_info[task_id].start = i;
		degree_sum_temp = 0;
		while((degree_sum_temp < max_degree_sum_task)&&(i < nb_node))
		{
			if(to_degree == 1)
				switch(ds)
				{
					//case S_LOCAL:	degree_sum_temp = degree_sum_temp + (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8)) get_con_size_sl(Nodes_sl, i);
					case S_LOCAL:	degree_sum_temp = degree_sum_temp + (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8));
							break;

//					case S_GLOBAL:	degree_sum_temp = degree_sum_temp + (unsigned long long)get_con_size_sg(Nodes_sg, i);
					case S_GLOBAL:	degree_sum_temp = degree_sum_temp + (unsigned long long)(Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8)) ;
							break;

					case S_NONE:	degree_sum_temp = degree_sum_temp + (unsigned long long) Nodes[i].con_size;
							break;
	
					default:
							printf("\n Unknown representation!");
				}
			else
				switch(ds)
				{
					//case S_LOCAL: degree_sum_temp = degree_sum_temp + (unsigned long long)get_from_size_sl(Nodes_sl, i);
					case S_LOCAL: degree_sum_temp = degree_sum_temp + (unsigned long long) (Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8));
						      break;

					//case S_GLOBAL:degree_sum_temp = degree_sum_temp + (unsigned long long)get_from_size_sg(Nodes_sg, i);
					case S_GLOBAL:degree_sum_temp = degree_sum_temp + (unsigned long long)(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8));
						      break;

					case S_NONE:  degree_sum_temp = degree_sum_temp + (unsigned long long) Nodes[i].from_size;
						      break;
					default:
							printf("\n Unknown representation!");
				}
			i++;
		}
		if(i >= (nb_node - 1))
			task_info[task_id].end = nb_node;
		else
			task_info[task_id].end = i+1;
		task_info[task_id].tid = task_id;
		task_info[task_id].degree_sum = degree_sum_temp;
		//printf("task %d begin at %d end at %d degree_sum = %lld\n",task_id, task_info[task_id].start, task_info[task_id].end, task_info[task_id].degree_sum);
		task_id++;
		i++;
		real_nb_task = task_id;
	}

	while(task_id < nb_task)
	{
		task_info[task_id].start = nb_node;
		task_info[task_id].end = nb_node;
		task_info[task_id].tid = task_id;
		task_info[task_id].degree_sum = 0;
		//printf("task %d begin at %d end at %d degree_sum = %lld\n",task_id, task_info[task_id].start, task_info[task_id].end, task_info[task_id].degree_sum);
		task_id++;
	}
		//printf("nb_task = %d\treal_nb_task = %d\tdegree_sum = %lld max_degree_sum_task = %lld\n",nb_task, real_nb_task,degree_sum, max_degree_sum_task);

	return task_info;
}

int Threads_Allocation_community_load(char* file_load_comm)
{

	// Stores thread's data		
	Threads_data = generate_community_task(file_load_comm);
        	
	return real_nb_task;
}

int Threads_Allocation_degree(int N, int scd, int to_degree_, int num_threads_, int chunk_, int ds)
{

	int i;
	
	if(scd == pthread_ || scd == pthread_2)
		Threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));

	// Stores thread's data		
	Threads_data = generate_task_index(to_degree_, num_threads_, chunk_, ds);	
        	
	return real_nb_task;
}
int Threads_Allocation_nodes(int N, int scd, int num_threads, int chunk)
{

	int i, nb_task = num_threads*chunk;
	double N_split =  (double) N / nb_task;
	
	// Allocate memory for threads
	if(scd == pthread_ || scd == pthread_2)
		Threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));

	// Stores thread's data		
	Threads_data = (Thread*)malloc(nb_task * sizeof(Thread));	
	
	// Split dataset into subsets, given to each thread
	Threads_data[0].tid = 0;
	Threads_data[0].start = 0;
	Threads_data[0].end = floor(N_split);

	for (i = 1; i < nb_task; i++)
	{
		Threads_data[i].tid = i;
		Threads_data[i].start = Threads_data[i - 1].end;
		if (i < (nb_task - 1))
		{
			Threads_data[i].end = Threads_data[i].start + floor(N_split);
		}
		else
		{
			Threads_data[i].end = N;
		}
	}
	
	return nb_task;
}

/***** Memory allocation - Initializations for Nodes *****/
 Node* Nodes_Allocation(int N)
{

	Node *Nodes;
	int i;
	Nodes = (Node*)malloc(N*sizeof(Node));
    
    for (i = 0; i < N; i++)
	{
		Nodes[i].con_size = 0;
		Nodes[i].from_size = 0;
		Nodes[i].From_id = (int*) malloc(sizeof(int));
        }	

    return Nodes;
}

/***** Read graph connections from txt file *****/	

void Read_from_txt_file(char* filename, Node* Nodes)
{
    
    FILE *fid;

    int from_idx, to_idx;
	int temp_size;

    //fid = fopen("web-Google.txt", "r");
    fid = fopen(filename, "r");
   	if (fid == NULL){printf("Error opening the file\n");}

	while (!feof(fid))
	{
		if (fscanf(fid,"%d\t%d\n", &from_idx, &to_idx))
		{
		//	printf("%d\t%d\n",from_idx,to_idx);
			Nodes[from_idx].con_size++;
			Nodes[to_idx].from_size++;
			temp_size = Nodes[to_idx].from_size;
			Nodes[to_idx].From_id = (int*) realloc(Nodes[to_idx].From_id, temp_size * sizeof(int));
			Nodes[to_idx].From_id[temp_size - 1] = from_idx; 
			degree_sum_to++;
			degree_sum_from++;
		}
	}
	//printf("End of connections insertion!\n");
	fclose(fid);
}

/***** Read P vector from txt file*****/	

void Read_P_from_txt_file()
{

	FILE *fid;
	double temp_P;
	int index = 0;

    fid = fopen("P.txt", "r");
   	if (fid == NULL){printf("Error opening the Probabilities file\n");}

	while (!feof(fid))
	{
		// P's values are double!
		if (fscanf(fid,"%lf\n", &temp_P))
		{
			Nodes[index].p_t1 = temp_P;
			index++;	   
		}
	}
	printf("End of P insertion!");

	fclose(fid);	
}


/***** Read E vector from txt file*****/	

void Read_E_from_txt_file()
{

	FILE *fid;
	double temp_E;
	int index = 0;
	
    fid = fopen("E.txt", "r");
   	if (fid == NULL){printf("Error opening the E file\n");}

	while (!feof(fid))
	{
		// E's values are double!
		if (fscanf(fid,"%lf\n", &temp_E))
		{
			Nodes[index].e = temp_E;
			index++;   
		}
	}
	printf("End of E insertion!");

	fclose(fid);	

}

/***** Create P and E with equal probability *****/

void Random_P_E()
{

   	int i;
    // Sum of P (it must be =1)
    double sum_P_1 = 0;
    // Sum of E (it must be =1)
    double sum_E_1 = 0; 
    
    
    // Arrays initialization
    for (i = 0; i < N; i++)
    {
        Nodes[i].p_t0 = 0;
        Nodes[i].p_t1 = 1;
        Nodes[i].p_t1 = (double) Nodes[i].p_t1 / N;

        sum_P_1 = sum_P_1 + Nodes[i].p_t1;
        
		Nodes[i].e = 1;
        Nodes[i].e = (double) Nodes[i].e / N;
        sum_E_1 = sum_E_1 + Nodes[i].e;
    }

    assert(sum_P_1 = 1);
    
    assert(sum_E_1 = 1);
}
    
void Random_P_E_sl(Node_succint_local* Nodes_sl, int N)
{

   	int i;
    // Sum of P (it must be =1)
    double sum_P_1 = 0;
    // Sum of E (it must be =1)
    double sum_E_1 = 0; 
    
    
    // Arrays initialization
    for (i = 0; i < N; i++)
    {
        Nodes_sl[i].p_t0 = 0;
        Nodes_sl[i].p_t1 = 1;
        Nodes_sl[i].p_t1 = (double) Nodes_sl[i].p_t1 / N;

        sum_P_1 = sum_P_1 + Nodes_sl[i].p_t1;
        
		Nodes_sl[i].e = 1;
        Nodes_sl[i].e = (double) Nodes_sl[i].e / N;
        sum_E_1 = sum_E_1 + Nodes_sl[i].e;
    }

    assert(sum_P_1 = 1);
    
    assert(sum_E_1 = 1);
}
    
void Random_P_E_sg(Node_succint_global* Nodes_sg, int N)
{

   	int i;
    // Sum of P (it must be =1)
    double sum_P_1 = 0;
    // Sum of E (it must be =1)
    double sum_E_1 = 0; 
    
    
    // Arrays initialization
    for (i = 0; i < N; i++)
    {
        Nodes_sg[i].p_t0 = 0;
        Nodes_sg[i].p_t1 = 1;
        Nodes_sg[i].p_t1 = (double) Nodes_sg[i].p_t1 / N;

        sum_P_1 = sum_P_1 + Nodes_sg[i].p_t1;
        
		Nodes_sg[i].e = 1;
        Nodes_sg[i].e = (double) Nodes_sg[i].e / N;
        sum_E_1 = sum_E_1 + Nodes_sg[i].e;
    }

    assert(sum_P_1 = 1);
    
    assert(sum_E_1 = 1);
}

void call_Local_Max_task(int nb_task)
{
	int i;
	for(i=0; i < nb_task; i++)
	{
#pragma omp task
		Local_Max_task(Threads_data[i]);
	}
}	

void call_Pagerank_Parallel_task(int nb_task)
{
	int i;
	for(i=0; i < nb_task; i++)
	{
#pragma omp task
		Pagerank_Parallel_task(Threads_data[i]);
	}
}	

void call_reinit_task(int nb_task)
{
	int i;
	for(i=0; i < nb_task; i++)
	{
#pragma omp task
		P_reinit_task(Threads_data[i]);
	}
}	

void* P_reinit_openmp_strategies(int end_index, int scd, Node* Nodes)
{
	/***************** for load informations issu ********************/
	int tid, i, nb_task = end_index, nb_node = end_index;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;               
	timeval_t t1, t2;
	timezone_t tz;
	switch(scd)
	{
		case taskloop_:
#pragma omp parallel
#pragma omp single
#pragma omp taskloop 
			for (int j = 0; j < nb_task; j++){
				P_reinit_task(Threads_data[j]);
			}
//	}
			break;
		case static_:
#pragma omp parallel for schedule(static,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       P_reinit_task(Threads_data[i]);
			break;
		case dynamic_:
#pragma omp parallel for schedule(dynamic,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       P_reinit_task(Threads_data[i]);
			break;
		case guided_:
#pragma omp parallel for schedule(guided,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       P_reinit_task(Threads_data[i]);
			break;
		case auto_:
#pragma omp parallel for schedule(auto) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       P_reinit_task(Threads_data[i]);
			break;
		case static_2:
#pragma omp parallel for schedule(static,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       P_reinit_task(Threads_data[i]);
			break;
		case dynamic_2:
#pragma omp parallel for schedule(dynamic,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       P_reinit_task(Threads_data[i]);
			break;
		case guided_2:
#pragma omp parallel for schedule(guided,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       P_reinit_task(Threads_data[i]);
			break;
		case auto_2:
#pragma omp parallel for schedule(auto) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       P_reinit_task(Threads_data[i]);
			break;
		case taskloop_2:
#pragma omp parallel
#pragma omp single
#pragma omp taskloop 
			for (int j = 0; j < nb_task; j++){
				P_reinit_task(Threads_data[j]);
			}
//	}
			break;
		default:
			printf("Unknown scheduling in init, verify it !!!!\n");
	}	
	return 0;
}

void* P_reinit_task(Thread thread_data)
{
	int tid, i, begin = thread_data.start, end = thread_data.end;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;               
	unsigned long long start_time, end_time;
	timeval_t t1, t2;
	timezone_t tz;

	tid = omp_get_thread_num();

	top_(&t1,&tz);
	for (i = begin; i < end; i++)
	{
		switch(DS)
		{	
			case S_NONE: Nodes[i].p_t0 = Nodes[i].p_t1;	
				     Nodes[i].p_t1 = 0;
				     break;
			case S_GLOBAL: Nodes_sg[i].p_t0 = Nodes_sg[i].p_t1;	
				       Nodes_sg[i].p_t1 = 0;
				      break;
			case S_LOCAL: Nodes_sl[i].p_t0 = Nodes_sl[i].p_t1;	
				      Nodes_sl[i].p_t1 = 0;
				      break;
			default: printf("Unknown graph representation");
		}

		if(to_degree == 1)
				switch(DS)
				{
					//case S_LOCAL:	degree = (unsigned long long) get_con_size_sl(Nodes_sl, i);
					case S_LOCAL:	degree = (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8));
							break;

					//case S_GLOBAL:	degree = (unsigned long long)get_con_size_sg(Nodes_sg, i);
					case S_GLOBAL:	degree = (unsigned long long)(Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8));
							break;

					case S_NONE:	degree = (unsigned long long) Nodes[i].con_size;
							break;
	
					default:
							printf("\n Unknown representation!");
				}
		else
				switch(DS)
				{
					//case S_LOCAL: degree = (unsigned long long) get_from_size_sl(Nodes_sl, i);
					case S_LOCAL: degree = (unsigned long long)(Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8));
						      break;

					//case S_GLOBAL:degree = (unsigned long long)  get_from_size_sg(Nodes_sg, i);
					case S_GLOBAL:degree = (unsigned long long)(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8));
						      break;

					case S_NONE:  degree = (unsigned long long) Nodes[i].from_size;
						      break;
					default:
							printf("\n Unknown representation!");
				}
		tab_load_reinit[tid].nb_node_sum += 1;
		tab_load_reinit[tid].degree_sum  += degree;
	}
	top_(&t2,&tz);
	temps = cpu_time_(t1,t2);
	//printf("%d\t%d\t%d\t%ld.%03ld\n",tid, i, degree, temps/1000,temps%1000);

	tab_load_reinit[tid].time        += temps;
	return 0;
}


void* P_reinit(void* arg)
{
	int tid, i;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;               
	unsigned long long start_time, end_time ;
	timeval_t t1, t2;
	timezone_t tz;

	Thread *thread_data = (Thread *)arg;
	tid = thread_data->tid;

	//start_time = _rdtsc () ;
	top_(&t1,&tz);
	for (i = thread_data->start; i < thread_data->end; i++)
	{
		switch(DS)
		{	
			case S_NONE: Nodes[i].p_t0 = Nodes[i].p_t1;	
				     Nodes[i].p_t1 = 0;
				     break;
			case S_GLOBAL: Nodes_sg[i].p_t0 = Nodes_sg[i].p_t1;	
				       Nodes_sg[i].p_t1 = 0;
				      break;
			case S_LOCAL: Nodes_sl[i].p_t0 = Nodes_sl[i].p_t1;	
				      Nodes_sl[i].p_t1 = 0;
				      break;
			default: printf("Unknown graph representation");
		}

		if(to_degree == 1)
				switch(DS)
				{
					//case S_LOCAL:	degree = (unsigned long long) get_con_size_sl(Nodes_sl, i);
					case S_LOCAL:	degree = (unsigned long long) (Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8));
							break;

					//case S_GLOBAL:	degree = (unsigned long long) get_con_size_sg(Nodes_sg, i);
					case S_GLOBAL:	degree = (unsigned long long)(Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8));
							break;

					case S_NONE:	degree = (unsigned long long) Nodes[i].con_size;
							break;
	
					default:
							printf("\n Unknown representation!");
				}
		else
				switch(DS)
				{
					//case S_LOCAL: degree = (unsigned long long) get_from_size_sl(Nodes_sl, i);
					case S_LOCAL: degree = (unsigned long long)(Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8));
						      break;

					//case S_GLOBAL:degree = (unsigned long long) get_from_size_sg(Nodes_sg, i);
					case S_GLOBAL:degree = (unsigned long long)(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8));
						      break;

					case S_NONE:  degree = (unsigned long long) Nodes[i].from_size;
						      break;
					default:
							printf("\n Unknown representation!");
				}
		tab_load_reinit[tid].nb_node_sum += 1;
		tab_load_reinit[tid].degree_sum  += degree;
	}
	top_(&t2,&tz);
	//end_time   = _rdtsc () ;
	temps = cpu_time_(t1,t2);
	//printf("%d\t%d\t%d\t%ld.%03ld\n",tid, i, degree, temps/1000,temps%1000);

	tab_load_reinit[tid].time        += temps;
	return 0;
}


void verif_nodes(Node* Nodes, int nb_node)
{
	int i, j, index;
	int tid = omp_get_thread_num();
	
	printf("\nTid = %d, Je vérifie les noeuds\n", tid);	
	for (i = 0; i < nb_node; i++) 
	{
		if (Nodes[i].from_size != 0)
		{
		    for (j = 0; j < Nodes[i].from_size; j++)
		    {
				index = Nodes[i].From_id[j];
				if(index > nb_node) 
					printf("tid = %d, i = %d, j = %d, index = %d, nb_node = %d \n",tid, i,j, index, nb_node);	
				//Nodes[i].p_t1 = Nodes[i].p_t1 + (double) Nodes[index].p_t0 / Nodes[index].con_size;
		    }
		}
	}
	printf("Tid = %d J'ai vérifié les noeuds\n", tid);	
}

/***** Main parallel algorithm *****/


void* Pagerank_Parallel_openmp_strategies(int end_index, int scd, Node* Nodes)
{
	/***************** for load informations issu ********************/
	int tid, nb_task = end_index, nb_node = end_index;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;               
	timeval_t t1, t2;
	timezone_t tz;
	/***************** for load informations issu ********************/

	int i, j, index;

	// Every thread will compute a local sum and add it
	// to the global one
	double temp_sum = 0;

	switch(scd)
	{
		case taskloop_:
#pragma omp parallel
#pragma omp single
#pragma omp taskloop 
			for (int j = 0; j < nb_task; j++)
			       Pagerank_Parallel_task(Threads_data[j]);
//	}
			break;
		case static_:
#pragma omp parallel for schedule(static,chunk) private(temps,tid,degree, index)
			for (i = 0; i < nb_task; i++)
			       Pagerank_Parallel_task(Threads_data[i]);
			break;
		case dynamic_:
#pragma omp parallel for schedule(dynamic,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Pagerank_Parallel_task(Threads_data[i]);
			break;
		case guided_:
#pragma omp parallel for schedule(guided,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Pagerank_Parallel_task(Threads_data[i]);
			break;
		case auto_:
#pragma omp parallel for schedule(auto) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Pagerank_Parallel_task(Threads_data[i]);
			break;
		case static_2:
#pragma omp parallel for schedule(static,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Pagerank_Parallel_task(Threads_data[i]);
			break;
		case dynamic_2:
#pragma omp parallel for schedule(dynamic,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Pagerank_Parallel_task(Threads_data[i]);
			break;
		case guided_2:
#pragma omp parallel for schedule(guided,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Pagerank_Parallel_task(Threads_data[i]);
			break;
		case auto_2:
#pragma omp parallel for schedule(auto) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Pagerank_Parallel_task(Threads_data[i]);
			break;
		case taskloop_2:
#pragma omp parallel
#pragma omp single
/*#pragma omp taskgroup
	{	
		int j;*/
#pragma omp taskloop 
			for (int j = 0; j < nb_task; j++)
			       Pagerank_Parallel_task(Threads_data[j]);
//	}
			break;

		default:
			printf("Unknown scheduling, verify it !!!!\n");

	}	
	return 0;
}

void* Pagerank_Parallel_task(Thread thread_data)
{
	/***************** for load informations issu ********************/
	int tid, begin = thread_data.start, end = thread_data.end;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;
	unsigned long long start_time, end_time;
	timeval_t t1, t2;
	timezone_t tz;
	/***************** for load informations issu ********************/

	int i, j, index;
	//tid = thread_data->tid;
	tid = omp_get_thread_num();

	// Every thread will compute a local sum and add it
	// to the global one
	double temp_sum = 0;
	int from_size;

	top_(&t1,&tz);
	for (i = begin; i < end; i++)
	{

		switch(DS)
		{
			case S_LOCAL:	
					//if((Nodes_sl[id_node].con_size[0] | (Nodes_sl[id_node].con_size[1]<<8))get_con_size_sl(Nodes_sl, i) == 0)
					if((Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8)) == 0)
					{
						 temp_sum = temp_sum + (double) Nodes_sl[i].p_t0 / ((double)N);
					}
					//if (get_from_size_sl(Nodes_sl, i) != 0)
					if ((Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8)) != 0)
					{
						from_size = (Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8)); 
					    // Compute the total probability, contributed by node's neighbors
					    for (j = 0; j < from_size; j++)
					    {
						//index = get_neighbor_sl(Nodes_sl, i, j);	
						if(Nodes_sl[i].neigh_size == 1)
							index = Nodes_sl[i].From_id[j][0];
						else 
							if(Nodes_sl[i].neigh_size == 2)
								index = (Nodes_sl[i].From_id[j][0] | (Nodes_sl[i].From_id[j][1]<<8));
							else
								index =  Nodes_sl[i].From_id[j][0]|(Nodes_sl[i].From_id[j][1]<<8)|(Nodes_sl[i].From_id[j][2]<<16);
						//Nodes_sl[i].p_t1 = Nodes_sl[i].p_t1 + (double) Nodes_sl[index].p_t0 /((double) get_con_size_sl(Nodes_sl, index));
						Nodes_sl[i].p_t1 = Nodes_sl[i].p_t1 + (double) Nodes_sl[index].p_t0 /((double) (Nodes_sl[index].con_size[0] | (Nodes_sl[index].con_size[1]<<8)));
					    }
					}
					break;

			case S_GLOBAL:	
					//if( get_con_size_sg(Nodes_sg, i) == 0)
					if((Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8)) == 0)
					{
						 temp_sum = temp_sum + (double) Nodes_sg[i].p_t0 / ((double)N);
					}

					//if (get_from_size_sg(Nodes_sg, i) != 0)
					if ((Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8)) != 0)
					{
						//from_size = get_from_size_sg(Nodes_sg, i); 
						from_size =(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8)); 
					    // Compute the total probability, contributed by node's neighbors
					    for (j = 0; j < from_size; j++)
					    {
						//index = get_neighbor_sg(Nodes_sg, i, j);	
						index = (Nodes_sg[i].From_id[j][0] | (Nodes_sg[i].From_id[j][1]<<8)| (Nodes_sg[i].From_id[j][2]<<16));	
						//Nodes_sg[i].p_t1 = Nodes_sg[i].p_t1 + (double) Nodes_sg[index].p_t0 /((double) get_con_size_sg(Nodes_sg, index));
						Nodes_sg[i].p_t1 = Nodes_sg[i].p_t1 + (double) Nodes_sg[index].p_t0 /((double)(Nodes_sg[index].con_size[0] | (Nodes_sg[index].con_size[1]<<8)));
					    }
					}
					break;

			case S_NONE:	
					if (Nodes[i].con_size == 0)
					{
						 temp_sum = temp_sum + (double) Nodes[i].p_t0 / ((double)N);
					}

					if (Nodes[i].from_size != 0)
					{
					    // Compute the total probability, contributed by node's neighbors
					    for (j = 0; j < Nodes[i].from_size; j++)
					    {
						index = Nodes[i].From_id[j];	
						Nodes[i].p_t1 = Nodes[i].p_t1 + (double) Nodes[index].p_t0 /((double) Nodes[index].con_size);
					    }
					}
					break;
	
			default:
					printf("\n Unknown representation!");
		}


		if(to_degree == 1)
				switch(DS)
				{
					//case S_LOCAL:	degree = (unsigned long long) get_con_size_sl(Nodes_sl, i);
					case S_LOCAL:	degree = (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8));
							break;

					//case S_GLOBAL:	degree = (unsigned long long) get_con_size_sg(Nodes_sg, i);
					case S_GLOBAL:	degree = (unsigned long long)(Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8));
							break;

					case S_NONE:	degree =  (unsigned long long) Nodes[i].con_size;
							break;
	
					default:
							printf("\n Unknown representation!");
				}
		else
				switch(DS)
				{
					//case S_LOCAL: degree = (unsigned long long) get_from_size_sl(Nodes_sl, i);
					case S_LOCAL: degree = (unsigned long long) (Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8));
						      break;

					//case S_GLOBAL:degree = (unsigned long long) get_from_size_sg(Nodes_sg, i);
					case S_GLOBAL:degree = (unsigned long long)(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8));
						      break;

					case S_NONE:  degree = (unsigned long long) Nodes[i].from_size;
						      break;
					default:
							printf("\n Unknown representation!");
				}
		tab_load_ppgr[tid].nb_node_sum += 1;
		tab_load_ppgr[tid].degree_sum  += degree;
	}
	
	// This is an atomic operation
	pthread_mutex_lock(&locksum);
	sum = sum + temp_sum; 
	pthread_mutex_unlock(&locksum);
/************************/		
	top_(&t2,&tz);
	temps = cpu_time_(t1,t2);
	tab_load_ppgr[tid].time        += temps;
	//printf("%d\t%d\t%d\t%ld.%03ld\n",tid, i, degree, temps/1000,temps%1000);
	return 0;
}


void* Pagerank_Parallel(void* arg)
{
	/***************** for load informations issu ********************/
	int tid;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;               
	unsigned long long start_time, end_time ;
	timeval_t t1, t2;
	timezone_t tz;
	/***************** for load informations issu ********************/

	Thread *thread_data = (Thread *) arg;
	int i, j, index, from_size;
	tid = thread_data->tid;

	// Every thread will compute a local sum and add it
	// to the global one
	double temp_sum = 0;

	top_(&t1,&tz);
	for (i = thread_data->start; i < thread_data->end; i++)
	{
		switch(DS)
		{
			case S_LOCAL:	
					//if(get_con_size_sl(Nodes_sl, i) == 0)
					if((Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8)) == 0)
					{
						 temp_sum = temp_sum + (double) Nodes_sl[i].p_t0 / ((double)N);
					}

					//if (get_from_size_sl(Nodes_sl, i) != 0)
					if ((Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8)) != 0)
					{
						from_size =(Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8)); 
					    // Compute the total probability, contributed by node's neighbors
					    for (j = 0; j < from_size; j++)
					    {
						//index = get_neighbor_sl(Nodes_sl, i, j);	
						if(Nodes_sl[i].neigh_size == 1)
							index = Nodes_sl[i].From_id[j][0];
						else 
							if(Nodes_sl[i].neigh_size == 2)
								index = (Nodes_sl[i].From_id[j][0] | (Nodes_sl[i].From_id[j][1]<<8));
							else
								index =  (Nodes_sl[i].From_id[j][0]|(Nodes_sl[i].From_id[j][1]<<8)|(Nodes_sl[i].From_id[j][2]<<16));
						
						//Nodes_sl[i].p_t1 = Nodes_sl[i].p_t1 + (double) Nodes_sl[index].p_t0 /((double) get_con_size_sl(Nodes_sl, index));
						Nodes_sl[i].p_t1 = Nodes_sl[i].p_t1 + (double) Nodes_sl[index].p_t0 /((double) (Nodes_sl[index].con_size[0] | (Nodes_sl[index].con_size[1]<<8)));
					    }
					}
					break;

			case S_GLOBAL:	
					//if(get_con_size_sg(Nodes_sg, i) == 0)
					if((Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8)) == 0)
					{
						 temp_sum = temp_sum + (double) Nodes_sg[i].p_t0 / ((double)N);
					}

					//if (get_from_size_sg(Nodes_sg, i) != 0)
					if ((Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8))!= 0)
					{
						//from_size = get_from_size_sg(Nodes_sg, i); 
						from_size =(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8)); 
					    // Compute the total probability, contributed by node's neighbors
					    for (j = 0; j < from_size; j++)
					    {
						//index =get_neighbor_sg(Nodes_sg, i, j);	
						index =(Nodes_sg[i].From_id[j][0] | (Nodes_sg[i].From_id[j][1]<<8)| (Nodes_sg[i].From_id[j][2]<<16));	
						//Nodes_sg[i].p_t1 = Nodes_sg[i].p_t1 + (double) Nodes_sg[index].p_t0 /((double) get_con_size_sg(Nodes_sg, index));
						Nodes_sg[i].p_t1 = Nodes_sg[i].p_t1 + (double) Nodes_sg[index].p_t0 /((double)(Nodes_sg[index].con_size[0] | (Nodes_sg[index].con_size[1]<<8)));
					    }
					}
					break;

			case S_NONE:	
					if (Nodes[i].con_size == 0)
					{
						 temp_sum = temp_sum + (double) Nodes[i].p_t0 / ((double)N);
					}

					if (Nodes[i].from_size != 0)
					{
					    // Compute the total probability, contributed by node's neighbors
					    for (j = 0; j < Nodes[i].from_size; j++)
					    {
						index = Nodes[i].From_id[j];	
						Nodes[i].p_t1 = Nodes[i].p_t1 + (double) Nodes[index].p_t0 /((double) Nodes[index].con_size);
					    }
					}
					break;
	
			default:
					printf("\n Unknown representation!");
		}


		if(to_degree == 1)
				switch(DS)
				{
					//case S_LOCAL:	degree = (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8)) get_con_size_sl(Nodes_sl, i);
					case S_LOCAL:	degree = (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8));
							break;

					//case S_GLOBAL:	degree = (unsigned long long) get_con_size_sg(Nodes_sg, i);
					case S_GLOBAL:	degree = (unsigned long long)(Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8));
							break;

					case S_NONE:	degree =  (unsigned long long) Nodes[i].con_size;
							break;
	
					default:
							printf("\n Unknown representation!");
				}
		else
				switch(DS)
				{
					//case S_LOCAL: degree = (unsigned long long) get_from_size_sl(Nodes_sl, i);
					case S_LOCAL: degree = (unsigned long long)(Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8));
						      break;

					//case S_GLOBAL:degree = (unsigned long long) get_from_size_sg(Nodes_sg, i);
					case S_GLOBAL:degree = (unsigned long long)(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8));
						      break;

					case S_NONE:  degree = (unsigned long long) Nodes[i].from_size;
						      break;
					default:
							printf("\n Unknown representation!");
				}

		tab_load_ppgr[tid].nb_node_sum += 1;
		tab_load_ppgr[tid].degree_sum  += degree;
	}
	
	// This is an atomic operation
	pthread_mutex_lock(&locksum);
	sum = sum + temp_sum; 
	pthread_mutex_unlock(&locksum);
/************************/		
	top_(&t2,&tz);
	temps = cpu_time_(t1,t2);
	tab_load_ppgr[tid].time        += temps;
	//printf("%d\t%d\t%d\t%ld.%03ld\n",tid, i, degree, temps/1000,temps%1000);
	return 0;
}

/***** Compute local max (thread's data max) *****/
void* Local_Max_task(Thread thread_data)
{
	/***************** for load informations issu ********************/
	int tid, begin = thread_data.start, end = thread_data.end;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;               
	unsigned long long start_time, end_time ;
	timeval_t t1, t2;
	timezone_t tz;
	/***************** for load informations issu ********************/

	int i, j;
	tid = omp_get_thread_num();
		
	// Every thread will find a local max and then check
	// if this is a global one
	double temp_max = -1;

	top_(&t1,&tz);
	for (i = begin; i < end; i++)
	{
		switch(DS)
		{
		case S_LOCAL:	Nodes_sl[i].p_t1 = d * (Nodes_sl[i].p_t1 + sum) + (1 - d) * Nodes_sl[i].e;
 
				if (fabs(Nodes_sl[i].p_t1 - Nodes_sl[i].p_t0) > temp_max)
				{
				    temp_max  = fabs(Nodes_sl[i].p_t1 - Nodes_sl[i].p_t0);
				}		
				break;

		case S_GLOBAL:  Nodes_sg[i].p_t1 = d * (Nodes_sg[i].p_t1 + sum) + (1 - d) * Nodes_sg[i].e;
 
				if (fabs(Nodes_sg[i].p_t1 - Nodes_sg[i].p_t0) > temp_max)
				{
				    temp_max  = fabs(Nodes_sg[i].p_t1 - Nodes_sg[i].p_t0);
				}		
				break;

		case S_NONE:	Nodes[i].p_t1 = d * (Nodes[i].p_t1 + sum) + (1 - d) * Nodes[i].e;
 
				if (fabs(Nodes[i].p_t1 - Nodes[i].p_t0) > temp_max)
				{
				    temp_max  = fabs(Nodes[i].p_t1 - Nodes[i].p_t0);
				}		
				break;
		default:
				printf("\n Unknown representation!");
		}
/************************/		
		if(to_degree == 1)
				switch(DS)
				{
					//case S_LOCAL:	degree = (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8)) get_con_size_sl(Nodes_sl, i);
					case S_LOCAL:	degree = (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8));
							break;

					//case S_GLOBAL:	degree = (unsigned long long)get_con_size_sg(Nodes_sg, i);
					case S_GLOBAL:	degree = (unsigned long long)(Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8));
							break;

					case S_NONE:	degree = (unsigned long long) Nodes[i].con_size;
							break;
	
					default:
							printf("\n Unknown representation!");
				}
		else
				switch(DS)
				{
					//case S_LOCAL: degree = (unsigned long long)(Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8)) get_from_size_sl(Nodes_sl, i);
					case S_LOCAL: degree = (unsigned long long)(Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8));
						      break;

					//case S_GLOBAL:degree = (unsigned long long) get_from_size_sg(Nodes_sg, i);
					case S_GLOBAL:degree = (unsigned long long)(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8));
						      break;

					case S_NONE:  degree = (unsigned long long) Nodes[i].from_size;
						      break;
					default:
							printf("\n Unknown representation!");
				}
		tab_load_lmax[tid].nb_node_sum += 1;
		tab_load_lmax[tid].degree_sum  += degree;
/************************/		
	}

	// Check if we have a new global max
	// This is an atomic operaiton
	pthread_mutex_lock(&lockmax);
	
	if (max_error  < temp_max)
	{			
		max_error = temp_max;		
	}	
	pthread_mutex_unlock(&lockmax);	
/************************/		
	top_(&t2,&tz);
	temps = cpu_time_(t1,t2);
	//printf("%d\t%d\t%d\t%ld.%03ld\n",tid, i, degree, temps/1000,temps%1000);

	tab_load_lmax[tid].time        += temps;
	return 0;
}

void* Local_Max_openmp_strategies(int end_index, int scd, Node* Nodes)
{
	/***************** for load informations issu ********************/
	int tid, nb_task = end_index, nb_node = end_index;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;               
	timeval_t t1, t2;
	timezone_t tz;
	/***************** for load informations issu ********************/

	int i, j;
		
	// Every thread will find a local max and then check
	// if this is a global one
	double temp_max = -1;
//printf("nb_node = %d\n",nb_node);
	switch(scd)
	{
		case taskloop_:
#pragma omp parallel
#pragma omp single
/*#pragma omp taskgroup
	{	
		int j;*/
#pragma omp taskloop 
			for (int j = 0; j < nb_task; j++)
			       Local_Max_task(Threads_data[j]);
//	}
			break;
		case static_:
#pragma omp parallel for schedule(static,chunk) private(temps,tid,degree, temp_max)
			for (i = 0; i < nb_task; i++)
			       Local_Max_task(Threads_data[i]);
			break;
		case dynamic_:
#pragma omp parallel for schedule(dynamic,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Local_Max_task(Threads_data[i]);
			break;
		case guided_:
#pragma omp parallel for schedule(guided,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Local_Max_task(Threads_data[i]);
			break;
		case auto_:
#pragma omp parallel for schedule(auto) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Local_Max_task(Threads_data[i]);
			break;
		case static_2:
#pragma omp parallel for schedule(static,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Local_Max_task(Threads_data[i]);
			break;
		case dynamic_2:
#pragma omp parallel for schedule(dynamic,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Local_Max_task(Threads_data[i]);
			break;
		case guided_2:
#pragma omp parallel for schedule(guided,chunk) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Local_Max_task(Threads_data[i]);
			break;
		case auto_2:
#pragma omp parallel for schedule(auto) private(temps,tid,degree)
			for (i = 0; i < nb_task; i++)
			       Local_Max_task(Threads_data[i]);
			break;
		case taskloop_2:
#pragma omp parallel
#pragma omp single
/*#pragma omp taskgroup
	{	
		int j;*/
#pragma omp taskloop 
			for (int j = 0; j < nb_task; j++)
			       Local_Max_task(Threads_data[j]);
//	}
			break;

		default:
			printf("Unknown scheduling, verify it !!!!\n");
	}	
	return 0;
}

void* Local_Max(void* arg)
{
	/***************** for load informations issu ********************/
	int tid;
	unsigned long time_up;
        unsigned long long time_sum = 0,degree_sum = 0, temps, nb_node_sum = 0, degree = 0, nb_node_tid;               
	unsigned long long start_time, end_time ;
	timeval_t t1, t2;
	timezone_t tz;
	/***************** for load informations issu ********************/

	Thread *thread_data = (Thread *) arg;
	int i, j;
	tid = thread_data->tid;
		
	// Every thread will find a local max and then check
	// if this is a global one
	double temp_max = -1;

	top_(&t1,&tz);
	for (i = thread_data->start; i < thread_data->end; i++)
	{
		switch(DS)
		{
		case S_LOCAL:	Nodes_sl[i].p_t1 = d * (Nodes_sl[i].p_t1 + sum) + (1 - d) * Nodes_sl[i].e;
 
				if (fabs(Nodes_sl[i].p_t1 - Nodes_sl[i].p_t0) > temp_max)
				{
				    temp_max  = fabs(Nodes_sl[i].p_t1 - Nodes_sl[i].p_t0);
				}		
				break;

		case S_GLOBAL:  Nodes_sg[i].p_t1 = d * (Nodes_sg[i].p_t1 + sum) + (1 - d) * Nodes_sg[i].e;
 
				if (fabs(Nodes_sg[i].p_t1 - Nodes_sg[i].p_t0) > temp_max)
				{
				    temp_max  = fabs(Nodes_sg[i].p_t1 - Nodes_sg[i].p_t0);
				}		
				break;

		case S_NONE:	Nodes[i].p_t1 = d * (Nodes[i].p_t1 + sum) + (1 - d) * Nodes[i].e;
 
				if (fabs(Nodes[i].p_t1 - Nodes[i].p_t0) > temp_max)
				{
				    temp_max  = fabs(Nodes[i].p_t1 - Nodes[i].p_t0);
				}		
				break;
		default:
				printf("\n Unknown representation!");
		}
/************************/		
		if(to_degree == 1)
				switch(DS)
				{
					//case S_LOCAL:	degree = (unsigned long long)get_con_size_sl(Nodes_sl, i);
					case S_LOCAL:	degree = (unsigned long long)(Nodes_sl[i].con_size[0] | (Nodes_sl[i].con_size[1]<<8));
							break;

					//case S_GLOBAL:	degree = (unsigned long long) get_con_size_sg(Nodes_sg, i);
					case S_GLOBAL:	degree = (unsigned long long)(Nodes_sg[i].con_size[0] | (Nodes_sg[i].con_size[1]<<8));
							break;

					case S_NONE:	degree = (unsigned long long) Nodes[i].con_size;
							break;
	
					default:
							printf("\n Unknown representation!");
				}
		else
				switch(DS)
				{
					//case S_LOCAL: degree = (unsigned long long) get_from_size_sl(Nodes_sl, i);
					case S_LOCAL: degree = (unsigned long long)(Nodes_sl[i].from_size[0] | (Nodes_sl[i].from_size[1]<<8));
						      break;

					//case S_GLOBAL:degree = (unsigned long long) get_from_size_sg(Nodes_sg, i);
					case S_GLOBAL:degree = (unsigned long long)(Nodes_sg[i].from_size[0] | (Nodes_sg[i].from_size[1]<<8));
						      break;

					case S_NONE:  degree = (unsigned long long) Nodes[i].from_size;
						      break;
					default:
							printf("\n Unknown representation!");
				}
		tab_load_lmax[tid].nb_node_sum += 1;
		tab_load_lmax[tid].degree_sum  += degree;
/************************/		
	}

	// Check if we have a new global max
	// This is an atomic operaiton
	pthread_mutex_lock(&lockmax);
	
	if (max_error  < temp_max)
	{			
		max_error = temp_max;		
	}	
	pthread_mutex_unlock(&lockmax);	
/************************/		
	top_(&t2,&tz);
	temps = cpu_time_(t1,t2);
	//printf("%d\t%d\t%d\t%ld.%03ld\n",tid, i, degree, temps/1000,temps%1000);

	tab_load_lmax[tid].time        += temps;
	return 0;
}

/***** Pagerank main algortihm *****/
void Pagerank_openmp_task(int nb_task)
{
 	/***** Start of algorithm *****/
	
    int i, j, index;
	omp_set_num_threads(num_threads);
	// Continue if we don't have convergence yet
    while (max_error > threshold)
    {
    	max_error = -1;
		sum = 0;

	// P array re-Initialization
#pragma omp parallel
	{
#pragma omp single
		call_reinit_task(nb_task);
	}
#pragma omp taskwait



        // Find P for each webpage
#pragma omp parallel
	{
#pragma omp single
		call_Pagerank_Parallel_task(nb_task);
	}
#pragma omp taskwait

		// Find local and global max
#pragma omp parallel
	{
#pragma omp single
		call_Local_Max_task(nb_task);
	}
#pragma omp taskwait
        
//        printf("Max Error in iteration %d = %4lf\n", iterations+1, max_error);
        iterations++;
    }
}

void Pagerank_openmp_strategies(int end_index, int scd, Node* Nodes)
{
 	/***** Start of algorithm *****/
	
    int i, j, index;
	omp_set_num_threads(num_threads);
	// Continue if we don't have convergence yet
    while (max_error > threshold)
    {
    	max_error = -1;
		sum = 0;

	// P array re-Initialization
	P_reinit_openmp_strategies(end_index, scd, Nodes);

        // Find P for each webpage
	Pagerank_Parallel_openmp_strategies(end_index,scd, Nodes);

		// Find local and global max
	Local_Max_openmp_strategies(end_index, scd, Nodes);

//        printf("Max Error in iteration %d = %4lf\n", iterations+1, max_error);
        iterations++;
    }
}

void Pagerank()
{
 	/***** Start of algorithm *****/
	
    int i, j, index;
	
	// Continue if we don't have convergence yet
    while (max_error > threshold)
    {
    	max_error = -1;
		sum = 0;

	// P array re-Initialization
        for (i = 0; i < num_threads; i++)
        {
			pthread_create(&Threads[i], NULL, &P_reinit,(void*) &Threads_data[i]);
	}

	// Wait for all threads to "catch" this point
	for (i = 0; i < num_threads; i++)
	{
		pthread_join(Threads[i], NULL);
        }


        // Find P for each webpage
        for (i = 0; i < num_threads; i++)
        {
            pthread_create(&Threads[i], NULL, &Pagerank_Parallel, (void*) &Threads_data[i]);   
        }

		for (i = 0; i < num_threads; i++)
		{
			pthread_join(Threads[i], NULL);
		}

		// Find local and global max
		for (i = 0; i < num_threads; i++)
        {
            pthread_create(&Threads[i], NULL, &Local_Max, (void*) &Threads_data[i]);   
        }

		for (i = 0; i < num_threads; i++)
		{
			pthread_join(Threads[i], NULL);
		}
        
//        printf("Max Error in iteration %d = %4lf\n", iterations+1, max_error);
        iterations++;
    }
}


/***** main function *****/   

int main(int argc, char** argv)
{

    struct timeval start, end;
    
    int i,j,k, end_index;
        unsigned long long time_sum = 0,degree_sum = 0, temps,init_time, work_time, nb_node_sum = 0, degree = 0, nb_node_tid;               
	unsigned long long start_time, end_time, temps1, temps2;
	timeval_t t1, t2;
	timezone_t tz;
	double totaltime;
	
	// Check input arguments
	if (argc < 10)
	{
		printf("# arguments:9 \n./pagerank file.txt N threshold=0.001 d=0.85 nb_thread scd chunk_size to_degree=0|1 file_load_comm(if scd = 0)\n");
		printf("\t\t\tNote that nb_tasks = nb_thread * chunk_size\n");
		return 0;
	}

	// get arguments 
	char* filename = (char*)malloc(sizeof(char)*100), *nom_fich = (char*)malloc(sizeof(char)*100), *file_load_comm = (char*)malloc(sizeof(char)*100);
	strcpy(filename, argv[1]);
	if(argc > 9)
		strcpy(file_load_comm, argv[9]);
	N = atoi(argv[2]);
	threshold = atof(argv[3]);
	d = atof(argv[4]); 
	num_threads = atoi(argv[5]);
	scd = atoi(argv[6]);
	chunk = atoi(argv[7]);
	to_degree = atoi(argv[8]);
	DS = atoi(argv[9]);

	// Check input arguments
	if ((num_threads < 1) || (num_threads > maxThreads)) 
	{
		printf("Threads number must be >= 1 and  <= %d!\n", maxThreads);
		exit(1);
	}
	// OR read probabilities from files
	//Read_P_from_txt_file();
	//Read_E_from_txt_file();

	switch(DS)
	{
		case S_LOCAL:	Nodes = Nodes_Allocation(N);
				Read_from_txt_file(filename, Nodes);

				//set_nodes_local(Nodes, Nodes_sl, N);
				Nodes_sl = init_nodes_local(Nodes, N);
				//printf("\nSize Nodes_sl=%lu, Nodes_sl[0]=%lu, \nNodes_sl[0].p_t0=%lu,\nNodes_sl[0].p_t1=%lu,\nNodes_sl[0].e=%lu, \nNodes_sl[0].neigh_size=%lu, \nNodes_sl[0].From_id=%lu, \nNodes_sl[0].con_size=%lu, \nNodes_sl[0].from_size=%lu", sizeof(Nodes_sl), sizeof(Nodes_sl[0]), sizeof(Nodes_sl[0].p_t0),sizeof(Nodes_sl[0].p_t1),sizeof(Nodes_sl[0].e),sizeof(Nodes_sl[0].neigh_size), sizeof(Nodes_sl[0].From_id), sizeof(Nodes_sl[0].con_size), sizeof(Nodes_sl[0].from_size));
				//printf("\n\nSize simple Nodes=%lu, Nodes[0]=%lu, \nNodes[0].p_t0=%lu,\nNodes[0].p_t1=%lu,\nNodes[0].e=%lu, \nNodes[0].From_id=%lu, \nNodes[0].con_size=%lu, \nNodes[0].from_size=%lu",sizeof(Nodes), sizeof(Nodes[0]), sizeof(Nodes[0].p_t0),sizeof(Nodes[0].p_t1),sizeof(Nodes[0].e),sizeof(Nodes[0].From_id), sizeof(Nodes[0].con_size), sizeof(Nodes[0].from_size));
				delete_nodes(Nodes, N);
				break;

		case S_GLOBAL:	Nodes = Nodes_Allocation(N);
				Read_from_txt_file(filename, Nodes);

				//set_nodes_global(Nodes, Nodes_sg, N);
				Nodes_sg = init_nodes_global(Nodes, N);
				//printf("\nSize Nodes_sg=%lu, Nodes_sg[0]=%lu, \nNodes_sg[0].p_t0=%lu,\nNodes_sg[0].p_t1=%lu,\nNodes_sg[0].e=%lu, \nNodes_sg[0].From_id=%lu, \nNodes_sg[0].con_size=%lu, \nNodes_sg[0].from_size=%lu", sizeof(Nodes_sg), sizeof(Nodes_sg[0]), sizeof(Nodes_sl[0].p_t0),sizeof(Nodes_sl[0].p_t1),sizeof(Nodes_sl[0].e),sizeof(Nodes_sg[0].From_id), sizeof(Nodes_sg[0].con_size), sizeof(Nodes_sg[0].from_size));
				//printf("\n\nSize simple Nodes=%lu, Nodes[0]=%lu, \nNodes[0].p_t0=%lu,\nNodes[0].p_t1=%lu,\nNodes[0].e=%lu, \nNodes[0].From_id=%lu, \nNodes[0].con_size=%lu, \nNodes[0].from_size=%lu",sizeof(Nodes), sizeof(Nodes[0]), sizeof(Nodes[0].p_t0),sizeof(Nodes[0].p_t1),sizeof(Nodes[0].e),sizeof(Nodes[0].From_id), sizeof(Nodes[0].con_size), sizeof(Nodes[0].from_size));

				delete_nodes(Nodes, N);
				break;

		case S_NONE:	Nodes = Nodes_Allocation(N);
				Read_from_txt_file(filename, Nodes);
				break;
	
		default:
				printf("\n Unknown representation!");

	}

				printf("\n First step is ok\n");
// /*
	print_max_degree(Nodes, N);
	start_time = _rdtsc () ;
	top_(&t1,&tz);
    switch(scd)
    {
	    case 0:        chunk = 1;
			   real_nb_task = Threads_Allocation_community_load(file_load_comm);
			   break;
	    case pthread_: chunk = 1;
			   real_nb_task = Threads_Allocation_nodes(N, scd, num_threads, chunk);
			   printf("Scheduling with node \t");
			   break;
	    case taskloop_:
	    case task_:
	    case static_ :
	    case dynamic_:
	    case guided_:
	    //case auto_:   real_nb_task = Threads_Allocation_nodes(N, scd, num_threads, 100*chunk);
	    case auto_:   real_nb_task = Threads_Allocation_nodes(N, scd, num_threads, chunk);
			  printf("Scheduling with node \t");
			  break;
	    case pthread_2: chunk = 1;
	    		  real_nb_task = Threads_Allocation_degree(N, scd, to_degree, num_threads, chunk, DS);
			  printf("Scheduling with degree \t");
			  break; 
	    case static_2:   
	    case dynamic_2:
	    case guided_2:
	    case auto_2:
	    case taskloop_2:
	    case task_2: real_nb_task = Threads_Allocation_degree(N, scd, to_degree, num_threads, chunk, DS);
			printf("Scheduling with degree \t");
			break; 
	default:
		printf("Unknown scheduling!! Choose between 1 and 12 \n");
		exit(-1);
    }
	end_index = real_nb_task;

	tab_load_reinit = init_tab_load_info(num_threads);
	tab_load_lmax = init_tab_load_info(num_threads);
	tab_load_ppgr = init_tab_load_info(num_threads);

	char* resl_split[100], *nom_fich_load_reinit = (char*) malloc(sizeof(char)*200),*nom_fich_load_lmax = (char*)malloc(sizeof(char)*200), *nom_fich_load_ppgr = (char*)malloc(sizeof(char)*200), *temp_load = (char*)malloc(sizeof(char)*200);

	size_t size = strlen(filename)-4;
	strncpy(temp_load,filename,size);
	sprintf(nom_fich_load_reinit,"reinit_load_file.txt",temp_load);
	sprintf(nom_fich_load_lmax,"lmax_load_file.txt",temp_load);
	sprintf(nom_fich_load_ppgr,"ppgr_load_file.txt",temp_load);
	//printf("%s\n%s\n%s\n%s\n%s\n%s\n",resl_split[0],filename,nom_fich_load_reinit,nom_fich_load_lmax,nom_fich_load_ppgr);

				printf("\n Second step is ok\n");
	switch(DS)
	{
		case S_NONE: Random_P_E();
			     break;
		case S_GLOBAL: Random_P_E_sg(Nodes_sg, N);
			      break;
		case S_LOCAL: Random_P_E_sl(Nodes_sl, N);
			      break;
		default: printf("Unknown graph representation");
	}

	top_(&t2,&tz);
	end_time   = _rdtsc () ;
	temps1 = end_time - start_time;//cpu_time_(t1,t2);
	init_time = cpu_time_(t1,t2);
	printf("%s\t%d\t%d\t%ld.%03ldms %lld\t",filename, real_nb_task,scd,init_time/1000,init_time%1000,temps1);


	start_time = _rdtsc ();
	top_(&t1,&tz);
    switch(scd)
    {
	    case 0: Pagerank_openmp_task(end_index);
		    break;
	    case pthread_: 
	    case pthread_2: Pagerank();
		   		break; 
	    case static_ :
	    case dynamic_:
	    case guided_:
	    case auto_:
	    case static_2:   
	    case dynamic_2:
	    case guided_2:
	    case taskloop_:
	    case taskloop_2:
	    case auto_2: Pagerank_openmp_strategies(end_index, scd, Nodes);
		         break;
	    case task_: 
	    case task_2: Pagerank_openmp_task(end_index);
		    break;
	    default:
		    printf("Unknown scheduling!! Choose between -1 and 13 \n");
		    exit(-1);
    }		  
	top_(&t2,&tz);
	end_time   = _rdtsc ();
	temps2 = end_time - start_time;//cpu_time_(t1,t2);
	work_time = cpu_time_(t1,t2);

    printf("%ld.%03ldms\t", init_time/1000,init_time%1000);
    printf("%ld.%03ldms\n", work_time/1000,work_time%1000);

        FILE* time_file = fopen("pagerank_execution_time.txt", "a");
        fprintf(time_file,"%s\t%d\t%d\t%d\t%d\t%d\t%d\t%ld.%03ld\t%ld.%03ld\t%lld\t%lld\n",filename,iterations,real_nb_task,num_threads,scd,chunk,DS,init_time/1000,init_time%1000,work_time/1000,work_time%1000,temps1,temps2);
	fclose(time_file);
   
	write_load_info(nom_fich_load_reinit, num_threads, tab_load_reinit, chunk, scd);
	write_load_info(nom_fich_load_lmax, num_threads, tab_load_lmax, chunk, scd);
	write_load_info(nom_fich_load_ppgr, num_threads, tab_load_ppgr, chunk, scd);
// * /
    return (EXIT_SUCCESS);
}
