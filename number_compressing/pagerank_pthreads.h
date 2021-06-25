#ifndef PAGERANK_PTHREADS 
#define PAGERANK_PTHREADS

#include <omp.h>
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

typedef struct Thread
{
	int tid;
	int start, end;
	unsigned long long degree_sum;
} Thread; 

/***** Struct used for Nodes data *****/
typedef unsigned char con_size_t[2];
typedef unsigned char from_size_t[2];
typedef unsigned char node_size_3_t[3];
typedef unsigned char* node_size_x_t;

typedef struct Node
{
	double p_t0;
	double p_t1;
	double e;
	int *From_id;
	int con_size;
	int from_size;
}Node;

#pragma pack(1)
typedef struct Node_succint_global
{
	double p_t0;
	double p_t1;
	double e;
	node_size_3_t *From_id;
	con_size_t con_size;
	from_size_t from_size;
}Node_succint_global;

typedef struct Node_succint_local
{
	double p_t0;
	double p_t1;
	double e;
	unsigned char neigh_size;
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

typedef thread_load_info_t* tab_load_info_t;
void* Pagerank_Parallel_openmp_strategies(int end_index, int scd, Node* Nodes);
void* Pagerank_Parallel_task(Thread thread_data);
void* Local_Max_task(Thread thread_data);
void* Local_Max_openmp_strategies(int end_index, int scd, Node* Nodes);
void* P_reinit_task(Thread thread_data);
void* P_reinit_openmp_strategies(int end_index, int scd, Node* Nodes);
void call_reinit_task(int nb_task);
void call_Pagerank_Parallel_task(int nb_task);
void call_Local_Max_task(int nb_task);
void Pagerank_openmp_strategies(int end_index, int scd, Node* Nodes);
void Pagerank_openmp_task(int nb_task);
void Pagerank();
void verif_nodes(Node* Nodes, int nb_node);

void print_max_degree(Node* Nodes, int N);

//begin succintness
void set_con_size_sl(Node_succint_local* Nodes_sl, int id_node, int con_size);
int get_con_size_sl(Node_succint_local* Nodes_sl, int id_node);
void set_from_size_sl(Node_succint_local* Nodes_sl, int id_node, int from_size);
int get_from_size_sl(Node_succint_local* Nodes_sl, int id_node);
int get_neighbor_sl(Node_succint_local* Nodes_sl, int id_node, int neigh_pos);
void allocate_from_id_array_sl(Node_succint_local* Nodes_sl, int id_node, int nb_neighbor, int max_neigh);
void set_neighbor_sl(Node_succint_local* Nodes_sl, int id_node, int neigh_pos, int neigh_id);
void set_neigh_size_sl(Node_succint_local* Nodes_sl, int id_node, int max_neigh);
int get_con_size_sg(Node_succint_global* Nodes_sg, int id_node);
void set_con_size_sg(Node_succint_global* Nodes_sg, int id_node, int con_size);
int get_from_size_sg(Node_succint_global* Nodes_sg, int id_node);
void set_from_size_sg(Node_succint_global* Nodes_sg, int id_node, int from_size);
int get_neighbor_sg(Node_succint_global* Nodes_sg, int id_node, int neigh_pos);
void set_neighbor_sg(Node_succint_global* Nodes_sg, int id_node, int neigh_pos, int neigh_id);
void allocate_from_id_array_sg(Node_succint_global* Nodes_sg, int id_node, int nb_neighbor);
int get_max_node(Node* Nodes, int id_node);
void set_nodes_local(Node* Nodes, Node_succint_local* Nodes_sl, int N);
void set_nodes_global(Node* Nodes, Node_succint_global* Nodes_sg, int N);
void delete_nodes(Node* Nodes, int N);
Node_succint_global* init_nodes_global(Node* Nodes, int N);
Node_succint_local* init_nodes_local(Node* Nodes, int N);
//end succintness
#endif
