#include"function.h"
#include "pagerank_pthreads.h"

void print_int(int*elt)
{
	printf("%d",*elt);
}

void print_string(char* elt)
{
	printf("%s",elt);
}

char *trim_space(char *str)
{
  char *end;

  // Trim leading space
  while(isspace(*str)) str++;

  if(*str == 0)  // All spaces?
    return str;

  // Trim trailing space
  end = str + strlen(str) - 1;
  while(end > str && isspace(*end)) end--;

  // Write new null terminator
  *(end+1) = 0;

  return str;
}

char** str_split(char* str, char* car, char*reslt[])
{
	char *str1, *str2, *token, *subtoken;
        char *saveptr1, *saveptr2;
        int j;
	
	for (j = 1, str1 = str; ; j++, str1 = NULL) {
               token = strtok_r(str1, car, &saveptr1);
               if (token == NULL)
                   break;
		reslt[j-1] = token;
          }
	
	return reslt;
}

//begin succintness
void set_con_size_sl(Node_succint_local* Nodes_sl, int id_node, int con_size)
{
	Nodes_sl[id_node].con_size[0] = 0;
	Nodes_sl[id_node].con_size[1] = 0;

	Nodes_sl[id_node].con_size[0] = Nodes_sl[id_node].con_size[0] | con_size;
	Nodes_sl[id_node].con_size[1] = Nodes_sl[id_node].con_size[1] | (con_size >> 8);

}

int get_con_size_sl(Node_succint_local* Nodes_sl, int id_node)
{
	/*int con_size = 0;
	con_size = con_size | Nodes_sl[id_node].con_size[0];
	con_size = con_size | (Nodes_sl[id_node].con_size[1]<<8);
	return con_size;*/

	return (Nodes_sl[id_node].con_size[0] | (Nodes_sl[id_node].con_size[1]<<8));
}

void set_from_size_sl(Node_succint_local* Nodes_sl, int id_node, int from_size)
{
	Nodes_sl[id_node].from_size[0] = 0;
	Nodes_sl[id_node].from_size[1] = 0;

	Nodes_sl[id_node].from_size[0] = Nodes_sl[id_node].from_size[0] | from_size;
	Nodes_sl[id_node].from_size[1] = Nodes_sl[id_node].from_size[1] | (from_size >> 8);
}

int get_from_size_sl(Node_succint_local* Nodes_sl, int id_node)
{
	/*int from_size = 0;
	from_size = from_size | Nodes_sl[id_node].from_size[0];
	from_size = from_size | (Nodes_sl[id_node].from_size[1]<<8);
	return from_size;*/

	return (Nodes_sl[id_node].from_size[0] | (Nodes_sl[id_node].from_size[1]<<8));
}

int get_neighbor_sl(Node_succint_local* Nodes_sl, int id_node, int neigh_pos)
{
	/*int neigh_id = 0;
	switch(Nodes_sl[id_node].neigh_size)
	{
		case 1:	neigh_id = neigh_id | Nodes_sl[id_node].From_id[neigh_pos][0];
			break;

		case 2:	neigh_id = neigh_id | Nodes_sl[id_node].From_id[neigh_pos][0];
			neigh_id = neigh_id | (Nodes_sl[id_node].From_id[neigh_pos][1]<<8);
			break;

		case 3:	neigh_id = neigh_id | Nodes_sl[id_node].From_id[neigh_pos][0];
			neigh_id = neigh_id | (Nodes_sl[id_node].From_id[neigh_pos][1]<<8);
			neigh_id = neigh_id | (Nodes_sl[id_node].From_id[neigh_pos][2]<<16);
			break;

		default:printf("Uknown size, %d, we considere neigh_size = 3 bytes\n",Nodes_sl[id_node].neigh_size);
			neigh_id = neigh_id | Nodes_sl[id_node].From_id[neigh_pos][0];
			neigh_id = neigh_id | (Nodes_sl[id_node].From_id[neigh_pos][1]<<8);
			neigh_id = neigh_id | (Nodes_sl[id_node].From_id[neigh_pos][2]<<16);

	}
	return neigh_id;*/
	switch(Nodes_sl[id_node].neigh_size)
	{
		case 1:	return Nodes_sl[id_node].From_id[neigh_pos][0];

		case 2:	return (Nodes_sl[id_node].From_id[neigh_pos][0] | (Nodes_sl[id_node].From_id[neigh_pos][1]<<8));

		case 3:	return Nodes_sl[id_node].From_id[neigh_pos][0]|(Nodes_sl[id_node].From_id[neigh_pos][1]<<8)|(Nodes_sl[id_node].From_id[neigh_pos][2]<<16);

		default:printf("Uknown size, %d, we considere neigh_size = 3 bytes\n",Nodes_sl[id_node].neigh_size);
			return Nodes_sl[id_node].From_id[neigh_pos][0]|(Nodes_sl[id_node].From_id[neigh_pos][1]<<8)|(Nodes_sl[id_node].From_id[neigh_pos][2]<<16);

	}
}

void allocate_from_id_array_sl(Node_succint_local* Nodes_sl, int id_node, int nb_neighbor, int max_neigh)
{
	int i;
	set_neigh_size_sl(Nodes_sl, id_node, max_neigh);
	Nodes_sl[id_node].From_id = (node_size_x_t *)malloc(sizeof(node_size_x_t)*nb_neighbor);
	for(i=0; i<nb_neighbor; i++)
		Nodes_sl[id_node].From_id[i] = (node_size_x_t)malloc(sizeof(unsigned char)*Nodes_sl[id_node].neigh_size);
}

void set_neighbor_sl(Node_succint_local* Nodes_sl, int id_node, int neigh_pos, int neigh_id)
{
	switch(Nodes_sl[id_node].neigh_size)
	{
		case 1:	Nodes_sl[id_node].From_id[neigh_pos][0] = 0;

			Nodes_sl[id_node].From_id[neigh_pos][0] = Nodes_sl[id_node].From_id[neigh_pos][0] | neigh_id;
			break;

		case 2:	Nodes_sl[id_node].From_id[neigh_pos][0] = 0;
	 		Nodes_sl[id_node].From_id[neigh_pos][1] = 0;

			Nodes_sl[id_node].From_id[neigh_pos][0] = Nodes_sl[id_node].From_id[neigh_pos][0] | neigh_id;
	 		Nodes_sl[id_node].From_id[neigh_pos][1] = Nodes_sl[id_node].From_id[neigh_pos][1] | (neigh_id>>8);
			break;

		case 3:	Nodes_sl[id_node].From_id[neigh_pos][0] = 0;
	 		Nodes_sl[id_node].From_id[neigh_pos][1] = 0;
	 		Nodes_sl[id_node].From_id[neigh_pos][2] = 0;

			Nodes_sl[id_node].From_id[neigh_pos][0] = Nodes_sl[id_node].From_id[neigh_pos][0] | neigh_id;
	 		Nodes_sl[id_node].From_id[neigh_pos][1] = Nodes_sl[id_node].From_id[neigh_pos][1] | (neigh_id>>8);
	 		Nodes_sl[id_node].From_id[neigh_pos][2] = Nodes_sl[id_node].From_id[neigh_pos][2] | (neigh_id>>16);
			break;

		default:printf("Uknown size, %d, we considere neigh_size = 3 bytes\n",Nodes_sl[id_node].neigh_size);
			Nodes_sl[id_node].From_id[neigh_pos][0] = 0;
	 		Nodes_sl[id_node].From_id[neigh_pos][1] = 0;
	 		Nodes_sl[id_node].From_id[neigh_pos][2] = 0;

			Nodes_sl[id_node].From_id[neigh_pos][0] = Nodes_sl[id_node].From_id[neigh_pos][0] | neigh_id;
	 		Nodes_sl[id_node].From_id[neigh_pos][1] = Nodes_sl[id_node].From_id[neigh_pos][1] | (neigh_id>>8);
	 		Nodes_sl[id_node].From_id[neigh_pos][2] = Nodes_sl[id_node].From_id[neigh_pos][2] | (neigh_id>>16);
	}
}

void set_neigh_size_sl(Node_succint_local* Nodes_sl, int id_node, int max_neigh)
{
	Nodes_sl[id_node].neigh_size = 1;
	if((max_neigh>>16)&0xFF)
		Nodes_sl[id_node].neigh_size = 3;
	else
		if((max_neigh>>8)&0xFF)
			Nodes_sl[id_node].neigh_size = 2;
}

int get_con_size_sg(Node_succint_global* Nodes_sg, int id_node)
{
	/*int con_size = 0;
	con_size = con_size | Nodes_sg[id_node].con_size[0];
	con_size = con_size | (Nodes_sg[id_node].con_size[1]<<8);
	return con_size;*/

	return Nodes_sg[id_node].con_size[0] | (Nodes_sg[id_node].con_size[1]<<8);
}

void set_con_size_sg(Node_succint_global* Nodes_sg, int id_node, int con_size)
{
	Nodes_sg[id_node].con_size[0] = 0;
	Nodes_sg[id_node].con_size[1] = 0;

	Nodes_sg[id_node].con_size[0] = Nodes_sg[id_node].con_size[0] | con_size;
	Nodes_sg[id_node].con_size[1] = Nodes_sg[id_node].con_size[1] | (con_size >> 8);
}



int get_from_size_sg(Node_succint_global* Nodes_sg, int id_node)
{
	/*int from_size = 0;
	from_size = from_size | Nodes_sg[id_node].from_size[0];
	from_size = from_size | (Nodes_sg[id_node].from_size[1]<<8);
	return from_size;*/

	return Nodes_sg[id_node].from_size[0] | (Nodes_sg[id_node].from_size[1]<<8);
}

void set_from_size_sg(Node_succint_global* Nodes_sg, int id_node, int from_size)
{
	Nodes_sg[id_node].from_size[0] = 0;
	Nodes_sg[id_node].from_size[1] = 0;

	Nodes_sg[id_node].from_size[0] = Nodes_sg[id_node].from_size[0] | from_size;
	Nodes_sg[id_node].from_size[1] = Nodes_sg[id_node].from_size[1] | (from_size >> 8);

}

int get_neighbor_sg(Node_succint_global* Nodes_sg, int id_node, int neigh_pos)
{
/*	int neigh_id = 0;
	neigh_id = neigh_id | Nodes_sg[id_node].From_id[neigh_pos][0];
	neigh_id = neigh_id | (Nodes_sg[id_node].From_id[neigh_pos][1]<<8);
	neigh_id = neigh_id | (Nodes_sg[id_node].From_id[neigh_pos][2]<<16);
	return neigh_id;*/
	return Nodes_sg[id_node].From_id[neigh_pos][0] | (Nodes_sg[id_node].From_id[neigh_pos][1]<<8)| (Nodes_sg[id_node].From_id[neigh_pos][2]<<16);
}

void set_neighbor_sg(Node_succint_global* Nodes_sg, int id_node, int neigh_pos, int neigh_id)
{
	 Nodes_sg[id_node].From_id[neigh_pos][0] = 0;
	 Nodes_sg[id_node].From_id[neigh_pos][1] = 0;
	 Nodes_sg[id_node].From_id[neigh_pos][2] = 0;

	
	 Nodes_sg[id_node].From_id[neigh_pos][0] = Nodes_sg[id_node].From_id[neigh_pos][0] | neigh_id;
	 Nodes_sg[id_node].From_id[neigh_pos][1] = Nodes_sg[id_node].From_id[neigh_pos][1] | (neigh_id>>8);
	 Nodes_sg[id_node].From_id[neigh_pos][2] = Nodes_sg[id_node].From_id[neigh_pos][2] | (neigh_id>>16);
}
	

void allocate_from_id_array_sg(Node_succint_global* Nodes_sg, int id_node, int nb_neighbor)
{
	Nodes_sg[id_node].From_id = (node_size_3_t *)malloc(sizeof(node_size_3_t)*nb_neighbor);
}

int get_max_node(Node* Nodes, int id_node)
{
	int i, max_node = Nodes[id_node].From_id[0];
	for(i=1; i<Nodes[id_node].from_size; i++)
	{
		if(max_node<Nodes[id_node].From_id[i])
			max_node = Nodes[id_node].From_id[i];
	}
	return max_node;
}

void print_max_degree(Node* Nodes, int N)
{
	int i, j, max_con_size = 0, max_from_size = 0;
	
	for(i=0;i<N;i++)
	{
		if(max_con_size < Nodes[i].con_size)
			max_con_size = Nodes[i].con_size;

		if(max_from_size < Nodes[i].from_size)
			max_from_size = Nodes[i].from_size;
	}

	printf("max_con_degree = %d\tmax_from_degree = %d\n",max_con_size, max_from_size);
}

Node_succint_local* init_nodes_local(Node* Nodes, int N)
{
	int i, j, max_neigh;
	Node_succint_local* Nodes_sl = (Node_succint_local*) malloc(sizeof(Node_succint_local)*N);

	for(i=0;i<N;i++)
	{
		Nodes_sl[i].p_t0 = Nodes[i].p_t0; 
		Nodes_sl[i].p_t1 = Nodes[i].p_t1; 
		Nodes_sl[i].e = Nodes[i].e;
		max_neigh = get_max_node(Nodes, i);
		set_con_size_sl(Nodes_sl, i, Nodes[i].con_size);
		set_from_size_sl(Nodes_sl, i, Nodes[i].from_size);

		allocate_from_id_array_sl(Nodes_sl, i, Nodes[i].from_size, max_neigh);
		for(j=0; j<Nodes[i].from_size; j++)
		{
			set_neighbor_sl(Nodes_sl, i, j, Nodes[i].From_id[j]);
		}
	}
	return Nodes_sl;
}

void set_nodes_local(Node* Nodes, Node_succint_local* Nodes_sl, int N)
{
	int i, j, max_neigh;
	Nodes_sl = (Node_succint_local*) malloc(sizeof(Node_succint_local)*N);

	for(i=0;i<N;i++)
	{
		Nodes_sl[i].p_t0 = Nodes[i].p_t0; 
		Nodes_sl[i].p_t1 = Nodes[i].p_t1; 
		Nodes_sl[i].e = Nodes[i].e;
		max_neigh = get_max_node(Nodes, i);
		set_con_size_sl(Nodes_sl, i, Nodes[i].con_size);
		set_from_size_sl(Nodes_sl, i, Nodes[i].from_size);

		allocate_from_id_array_sl(Nodes_sl, i, Nodes[i].from_size, max_neigh);
		for(j=0; j<Nodes[i].from_size; j++)
		{
			set_neighbor_sl(Nodes_sl, i, j, Nodes[i].From_id[j]);
		}
	}
}

Node_succint_global* init_nodes_global(Node* Nodes, int N)
{
	int i, j, max_neigh;
	Node_succint_global* Nodes_sg = (Node_succint_global*) malloc(sizeof(Node_succint_global)*N);

	for(i=0;i<N;i++)
	{
		Nodes_sg[i].p_t0 = Nodes[i].p_t0; 
		Nodes_sg[i].p_t1 = Nodes[i].p_t1; 
		Nodes_sg[i].e = Nodes[i].e;
		set_con_size_sg(Nodes_sg, i, Nodes[i].con_size);
		set_from_size_sg(Nodes_sg, i, Nodes[i].from_size);

		allocate_from_id_array_sg(Nodes_sg, i, Nodes[i].from_size);
		for(j=0; j<Nodes[i].from_size; j++)
		{
			set_neighbor_sg(Nodes_sg, i, j, Nodes[i].From_id[j]);
		}
	}
	return Nodes_sg;
}

void set_nodes_global(Node* Nodes, Node_succint_global* Nodes_sg, int N)
{
	int i, j, max_neigh;
	Nodes_sg = (Node_succint_global*) malloc(sizeof(Node_succint_global)*N);

	for(i=0;i<N;i++)
	{
		Nodes_sg[i].p_t0 = Nodes[i].p_t0; 
		Nodes_sg[i].p_t1 = Nodes[i].p_t1; 
		Nodes_sg[i].e = Nodes[i].e;
		set_con_size_sg(Nodes_sg, i, Nodes[i].con_size);
		set_from_size_sg(Nodes_sg, i, Nodes[i].from_size);

		allocate_from_id_array_sg(Nodes_sg, i, Nodes[i].from_size);
		for(j=0; j<Nodes[i].from_size; j++)
		{
			set_neighbor_sg(Nodes_sg, i, j, Nodes[i].From_id[j]);
		}
	}
}

void delete_nodes(Node* Nodes, int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
		free(Nodes[i].From_id);
    }
    free(Nodes);
}
//end succintness
