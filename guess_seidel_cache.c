// ./pagerank web-Google.txt 916428 0.00001 0.85 8

/********************************************************************/
/*    Pagerank project 2014 - Parallel version                      */
/*    	*based on Cleve Moler's matlab implementation               */
/*                                                                  */
/*    Implemented by Nikos Katirtzis (nikos912000)                  */
/********************************************************************/

/******************** Includes - Defines ****************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdint.h>

// Cache line size (64 bytes on most modern processors)
#define CACHE_LINE_SIZE 64

/***** Struct for timestamps *****/
struct timeval start,end;

// Thread-local storage for intermediate results
typedef struct {
    atomic_uintptr_t local_sum;
    atomic_uintptr_t local_max_error;
    atomic_int local_converged;
    char padding[32];  // Fixed padding size
} ThreadLocal;

/***** Struct used for Threads data *****/
typedef struct {
    int tid;
    int start, end;
    ThreadLocal local;
    char padding[32];  // Fixed padding size
} Thread;

/***** Struct used for Nodes data *****/
typedef struct {
    double p_t0;
    double p_t1;
    double e;
    int con_size;
    int from_size;
    int *From_id;
    char padding[32];  // Fixed padding size
} Node;

#define maxThreads 64

/******************** Defines ****************/
// Number of nodes
int N, num_threads;

// Convergence threashold and algorithm's parameter d  
double threshold, d;

//Table of threads
pthread_t *Threads;

// Table with thread's data
Thread *Threads_data;

// Table of node's data
Node *Nodes;

// Add atomic operations for shared variables
atomic_int global_converged = 0;
atomic_uintptr_t global_max_error = 0;
atomic_uintptr_t global_sum = 0;

// Helper functions for atomic double operations
static inline double atomic_load_double(atomic_uintptr_t *ptr) {
    uintptr_t val = atomic_load(ptr);
    return *(double*)&val;
}

static inline void atomic_store_double(atomic_uintptr_t *ptr, double value) {
    uintptr_t val = *(uintptr_t*)&value;
    atomic_store(ptr, val);
}

/***** Memory allocation - Initializations for Threads *****/

void Threads_Allocation()
{
    int i;
    double N_split = (double)N / num_threads;
    
    // Allocate memory for threads with proper alignment
    Threads = (pthread_t *)aligned_alloc(CACHE_LINE_SIZE, num_threads * sizeof(pthread_t));
    Threads_data = (Thread *)aligned_alloc(CACHE_LINE_SIZE, num_threads * sizeof(Thread));
    
    // Split dataset into subsets, given to each thread
    Threads_data[0].tid = 0;
    Threads_data[0].start = 0;
    Threads_data[0].end = floor(N_split);
    Threads_data[0].local.local_sum = 0;
    Threads_data[0].local.local_max_error = 0;
    Threads_data[0].local.local_converged = 1;

    for (i = 1; i < num_threads; i++) {
        Threads_data[i].tid = i;
        Threads_data[i].start = Threads_data[i-1].end;
        Threads_data[i].local.local_sum = 0;
        Threads_data[i].local.local_max_error = 0;
        Threads_data[i].local.local_converged = 1;
        
        if (i < (num_threads - 1)) {
            Threads_data[i].end = Threads_data[i].start + floor(N_split);
        } else {
            Threads_data[i].end = N;
        }
    }
    
    printf("\n");

    for (i = 0; i < num_threads; i++)
    {
        printf("Thread %d, start = %d, end = %d\n", Threads_data[i].tid, Threads_data[i].start, Threads_data[i].end);
    }

    printf("\n");
}

/***** Memory allocation - Initializations for Nodes *****/
void Nodes_Allocation()
{

	int i;
	Nodes = (Node*)malloc(N*sizeof(Node));
    
    for (i = 0; i < N; i++)
	{
		Nodes[i].con_size = 0;
		Nodes[i].from_size = 0;
        Nodes[i].From_id = (int*) malloc(sizeof(int));
    }	

}

/***** Read graph connections from txt file *****/	

void Read_from_txt_file(char* filename)
{
    FILE *fid;
    int from_idx, to_idx;
    int temp_size;
    char line[1000];

    fid = fopen("web-Google.txt", "r");
    if (fid == NULL) {
        printf("Error opening the file\n");
        return;
    }

    while (fgets(line, sizeof(line), fid) != NULL) {
        // ignore sentences starting from #
        if (sscanf(line,"%d\t%d\n", &from_idx, &to_idx)) {
            Nodes[from_idx].con_size++;
            Nodes[to_idx].from_size++;
            temp_size = Nodes[to_idx].from_size;
            Nodes[to_idx].From_id = (int*) realloc(Nodes[to_idx].From_id, temp_size * sizeof(int));
            Nodes[to_idx].From_id[temp_size - 1] = from_idx; 
        }
    }

    printf("End of connections insertion!\n");
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
    double sum_P_1 = 0;
    double sum_E_1 = 0; 
    
    // Arrays initialization
    for (i = 0; i < N; i++) {
        Nodes[i].p_t0 = 1.0 / N;
        Nodes[i].p_t1 = 1.0 / N;
        sum_P_1 += Nodes[i].p_t1;
        
        Nodes[i].e = 1.0 / N;
        sum_E_1 += Nodes[i].e;
    }

    // Assert sum of probabilities is =1
    assert(fabs(sum_P_1 - 1.0) < 1e-10);
    assert(fabs(sum_E_1 - 1.0) < 1e-10);
}

/***** Re-initialize P(t) and P(t + 1) values *****/

void* P_reinit(void* arg)
{

	Thread *thread_data = (Thread *)arg;
	int i;

	for (i = thread_data->start; i < thread_data->end; i++)
	{
			Nodes[i].p_t0 = Nodes[i].p_t1;	
			Nodes[i].p_t1 = 0;
	}
	return 0;
}

/***** Main parallel algorithm *****/

void* Pagerank_Parallel(void* arg)
{
    Thread *thread_data = (Thread *) arg;
    int i, j, index;
    const int BLOCK_SIZE = 64; // Cache block size
    double local_sum = 0;
    double local_max_error = 0;
    int local_converged = 1;

    while (!atomic_load(&global_converged)) {
        // Reset local values
        local_sum = 0;
        local_max_error = 0;
        local_converged = 1;

        // Process local nodes in blocks
        for (i = thread_data->start; i < thread_data->end; i += BLOCK_SIZE) {
            int block_end = (i + BLOCK_SIZE < thread_data->end) ? i + BLOCK_SIZE : thread_data->end;
            
            // Prefetch next block
            if (i + BLOCK_SIZE < thread_data->end) {
                __builtin_prefetch(&Nodes[i + BLOCK_SIZE], 0, 0);
            }

            // Process current block
            for (int k = i; k < block_end; k++) {
                double old_p = Nodes[k].p_t0;
                double new_p = 0;

                // Handle dangling nodes
                if (Nodes[k].con_size == 0) {
                    local_sum += old_p / N;
                }

                // Compute new PageRank value with prefetching
                if (Nodes[k].from_size != 0) {
                    // Prefetch next node's data
                    if (k + 1 < block_end) {
                        __builtin_prefetch(&Nodes[k + 1], 0, 0);
                    }

                    for (j = 0; j < Nodes[k].from_size; j++) {
                        index = Nodes[k].From_id[j];
                        // Prefetch next From_id element
                        if (j + 1 < Nodes[k].from_size) {
                            __builtin_prefetch(&Nodes[k].From_id[j + 1], 0, 0);
                        }
                        new_p += Nodes[index].p_t0 / Nodes[index].con_size;
                    }
                }

                // Apply damping factor and teleportation
                new_p = d * (new_p + atomic_load_double(&global_sum)) + (1 - d) * Nodes[k].e;

                // Update node's PageRank
                Nodes[k].p_t0 = new_p;
                Nodes[k].p_t1 = new_p;

                // Check convergence
                double error = fabs(new_p - old_p);
                if (error > local_max_error) {
                    local_max_error = error;
                }
                if (error > threshold) {
                    local_converged = 0;
                }
            }
        }

        // Update thread-local values atomically
        atomic_store_double(&thread_data->local.local_sum, local_sum);
        atomic_store_double(&thread_data->local.local_max_error, local_max_error);
        atomic_store(&thread_data->local.local_converged, local_converged);

        // Check global convergence
        int all_converged = 1;
        for (i = 0; i < num_threads; i++) {
            if (!atomic_load(&Threads_data[i].local.local_converged)) {
                all_converged = 0;
                break;
            }
        }

        if (all_converged) {
            atomic_store(&global_converged, 1);
        }

        // Update global sum
        double total_sum = 0;
        for (i = 0; i < num_threads; i++) {
            total_sum += atomic_load_double(&Threads_data[i].local.local_sum);
        }
        atomic_store_double(&global_sum, total_sum);

        // Update global max error
        double max_error = 0;
        for (i = 0; i < num_threads; i++) {
            double thread_error = atomic_load_double(&Threads_data[i].local.local_max_error);
            if (thread_error > max_error) {
                max_error = thread_error;
            }
        }
        atomic_store_double(&global_max_error, max_error);
    }

    return NULL;
}

/***** Compute local max (thread's data max) *****/
void* Local_Max(void* arg)
{
    Thread *thread_data = (Thread *) arg;
    int i;
    double local_max = 0;

    for (i = thread_data->start; i < thread_data->end; i++) {
        double error = fabs(Nodes[i].p_t1 - Nodes[i].p_t0);
        if (error > local_max) {
            local_max = error;
        }
    }

    atomic_store_double(&thread_data->local.local_max_error, local_max);
    return NULL;
}

/***** Pagerank main algortihm *****/
void Pagerank()
{
    int i;
    
    // Initialize global state
    atomic_store(&global_converged, 0);
    atomic_store_double(&global_max_error, 1.0);
    atomic_store_double(&global_sum, 0.0);

    // Create threads
    for (i = 0; i < num_threads; i++) {
        pthread_create(&Threads[i], NULL, &Pagerank_Parallel, (void*) &Threads_data[i]);
    }

    // Wait for all threads to complete
    for (i = 0; i < num_threads; i++) {
        pthread_join(Threads[i], NULL);
    }

    printf("Final max error: %f\n", atomic_load_double(&global_max_error));
}


/***** main function *****/   

int main(int argc, char** argv)
{
    struct timeval start, end;
    int i;
    double totaltime;
    
    // Check input arguments
    if (argc < 5)
    {
        printf("Error in arguments! Three arguments required: graph filename, N, threshold and d\n");
        return 0;
    }

    // get arguments 
    char filename[256];
    strcpy(filename, argv[1]);
    N = atoi(argv[2]);
    threshold = atof(argv[3]);
    d = atof(argv[4]); 
    num_threads = atoi(argv[5]);

    // Check input arguments
    if ((num_threads < 1) || (num_threads > maxThreads)) 
    {
        printf("Threads number must be >= 1 and  <= %d!\n", maxThreads);
        exit(1);
    }

    Threads_Allocation();
    Nodes_Allocation();
    
    // OR read probabilities from files
    Read_from_txt_file(filename);
    //Read_P_from_txt_file();
    //Read_E_from_txt_file();

    Random_P_E();

    printf("\n");
    printf("Parallel version of Pagerank\n");

    gettimeofday(&start, NULL);
    Pagerank();
    gettimeofday(&end, NULL);  

    // Print top 10 pages by rank
    printf("\nTop 10 pages by rank:\n");
    
    // Create array of indices for sorting
    int *indices = (int*)malloc(N * sizeof(int));
    for(i = 0; i < N; i++) {
        indices[i] = i;
    }
    
    // Sort indices based on page rank values
    for(i = 0; i < 10; i++) {
        for(int j = i+1; j < N; j++) {
            if(Nodes[indices[j]].p_t1 > Nodes[indices[i]].p_t1) {
                int temp = indices[i];
                indices[i] = indices[j]; 
                indices[j] = temp;
            }
        }
        printf("%d. Page %d (rank: %f)\n", i+1, indices[i], Nodes[indices[i]].p_t1);
    }
    
    free(indices);

    totaltime = (((end.tv_usec - start.tv_usec) / 1.0e6 + end.tv_sec - start.tv_sec) * 1000) / 1000;

    printf("\nTotaltime = %f seconds\n", totaltime);
    printf("End of program!\n");
    
    return (EXIT_SUCCESS);
}
