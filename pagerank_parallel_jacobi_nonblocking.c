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
#include <unistd.h>  // for usleep
#include <stdatomic.h>


/***** Struct for timestamps *****/
struct timeval start,end;

/***** Struct used for Threads data *****/

typedef struct
{
	int tid;
	int start, end;
} Thread; 

/***** Struct used for Nodes data *****/

typedef struct
{
	double p_t0;
	double p_t1;
	double e;
	int *From_id;
	int con_size;
	int from_size;
}Node;

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
pthread_mutex_t lockP = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t locksum = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t lockmax = PTHREAD_MUTEX_INITIALIZER;

// Number of threads
int num_threads;

// Number of iterations
int iterations = 0;

// Update shared variables to be atomic
_Atomic double sum = 0;
_Atomic double max_error = 1;
_Atomic int active_threads = 0;

// Add thread completion tracking
_Atomic int completed_threads = 0;

/***** Memory allocation - Initializations for Threads *****/

void Threads_Allocation()
{

	int i;
	double N_split =  (double) N / num_threads;
	
	// Allocate memory for threads
	Threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));

	// Stores thread's data		
	Threads_data = (Thread*)malloc(num_threads * sizeof(Thread));	
	
	// Split dataset into subsets, given to each thread
	Threads_data[0].tid = 0;
	Threads_data[0].start = 0;
	Threads_data[0].end = floor(N_split);

	for (i = 1; i < num_threads; i++)
	{
		Threads_data[i].tid = i;
		Threads_data[i].start = Threads_data[i - 1].end;
		if (i < (num_threads - 1))
		{
			Threads_data[i].end = Threads_data[i].start + floor(N_split);
		}
		else
		{
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
    if (fid == NULL){printf("Error opening the file\n");}

    while (!feof(fid))
    {
        if (fgets(line, sizeof(line), fid) == NULL) {
            break;
        }
        // ignore sentences starting from #
        if (sscanf(line,"%d\t%d\n", &from_idx, &to_idx))
        {
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

    // Assert sum of probabilities is =1
    
    // Print sum of P (it must be =1)
    //printf("Sum of P = %f\n",sum_P_1);
    
    // Exit if sum of P is !=1
    assert(sum_P_1 = 1);
    
    //printf("\n");
    
    // Print sum of E (it must be =1)
    //printf("Sum of E = %f\n",sum_E_1);
    
    // Exit if sum of Pt0 is !=1
    assert(sum_E_1 = 1);

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
    atomic_fetch_add(&completed_threads, 1);
    return 0;
}

/***** Main parallel algorithm *****/

void* Pagerank_Parallel(void* arg)
{
    Thread *thread_data = (Thread *) arg;
    int i, j, index;
    double temp_sum = 0;

    for (i = thread_data->start; i < thread_data->end; i++)
    {
        if (Nodes[i].con_size == 0)
        {
            temp_sum = temp_sum + (double) Nodes[i].p_t0 / N;
        }

        if (Nodes[i].from_size != 0)
        {
            for (j = 0; j < Nodes[i].from_size; j++)
            {
                index = Nodes[i].From_id[j];    
                Nodes[i].p_t1 = Nodes[i].p_t1 + (double) Nodes[index].p_t0 / Nodes[index].con_size;
            }
        }        
    }
    
    // Atomic add to global sum using atomic_fetch_add_explicit
    double current_sum;
    do {
        current_sum = atomic_load_explicit(&sum, memory_order_relaxed);
    } while (!atomic_compare_exchange_weak_explicit(&sum, &current_sum, current_sum + temp_sum,
                                                  memory_order_release, memory_order_relaxed));
    atomic_fetch_add(&completed_threads, 1);
    return 0;
}

/***** Compute local max (thread's data max) *****/
void* Local_Max(void* arg)
{
    Thread *thread_data = (Thread *) arg;
    int i, j;
    double temp_max = -1;

    for (i = thread_data->start; i < thread_data->end; i++)
    {
        double current_sum = atomic_load_explicit(&sum, memory_order_acquire);
        Nodes[i].p_t1 = d * (Nodes[i].p_t1 + current_sum) + (1 - d) * Nodes[i].e;
 
        if (fabs(Nodes[i].p_t1 - Nodes[i].p_t0) > temp_max)
        {
            temp_max = fabs(Nodes[i].p_t1 - Nodes[i].p_t0);
        }
    }

    // Atomic compare and swap for max_error
    double current_max;
    do {
        current_max = atomic_load_explicit(&max_error, memory_order_relaxed);
        if (temp_max <= current_max) break;
    } while (!atomic_compare_exchange_weak_explicit(&max_error, &current_max, temp_max,
                                                  memory_order_release, memory_order_relaxed));
    
    atomic_fetch_add(&completed_threads, 1);
    return 0;
}

/***** Pagerank main algortihm *****/
void Pagerank()
{
    int i;
    pthread_t *threads_phase1 = malloc(num_threads * sizeof(pthread_t));
    pthread_t *threads_phase2 = malloc(num_threads * sizeof(pthread_t));
    pthread_t *threads_phase3 = malloc(num_threads * sizeof(pthread_t));
    
    while (atomic_load_explicit(&max_error, memory_order_acquire) > threshold)
    {
        atomic_store_explicit(&max_error, -1, memory_order_release);
        atomic_store_explicit(&sum, 0, memory_order_release);
        atomic_store(&completed_threads, 0);
        
        // Phase 1: Reinitialization
        for (i = 0; i < num_threads; i++)
        {
            pthread_create(&threads_phase1[i], NULL, &P_reinit, (void*) &Threads_data[i]);
        }
        
        // Wait for phase 1 to complete
        while (atomic_load(&completed_threads) < num_threads) {
            usleep(25);
        }
        
        atomic_store(&completed_threads, 0);
        
        // Phase 2: PageRank computation
        for (i = 0; i < num_threads; i++)
        {
            pthread_create(&threads_phase2[i], NULL, &Pagerank_Parallel, (void*) &Threads_data[i]);   
        }
        
        // Wait for phase 2 to complete
        while (atomic_load(&completed_threads) < num_threads) {
            usleep(25);
        }
        
        atomic_store(&completed_threads, 0);
        
        // Phase 3: Local max computation
        for (i = 0; i < num_threads; i++)
        {
            pthread_create(&threads_phase3[i], NULL, &Local_Max, (void*) &Threads_data[i]);   
        }
        
        // Wait for phase 3 to complete
        while (atomic_load(&completed_threads) < num_threads) {
            usleep(25);
        }
        
        printf("Max Error in iteration %d = %f\n", iterations+1, atomic_load_explicit(&max_error, memory_order_acquire));
        iterations++;
    }
    
    free(threads_phase1);
    free(threads_phase2);
    free(threads_phase3);
}


/***** main function *****/   

int main(int argc, char** argv)
{

    struct timeval start, end;
    
    int i,j,k;
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

    /*for (i = 0; i < N; i++)
    {
        printf("P_t1[%d] = %f\n",i, Nodes[i].p_t1);
    }
    
    printf("\n");*/
    
	// Print no of iterations
    printf("Total iterations: %d\n", iterations);


    totaltime = (((end.tv_usec - start.tv_usec) / 1.0e6+ end.tv_sec - start.tv_sec) * 1000) / 1000;

	printf("\nTotaltime = %f seconds\n", totaltime);
    
    printf("End of program!\n");
    
    return (EXIT_SUCCESS);
}
