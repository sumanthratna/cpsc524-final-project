/******************** Includes - Defines ****************/
#include <assert.h>
#include <math.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>
#include <limits.h>

/***** Struct for timestamps *****/
struct timeval start, end;

/***** Struct used for Threads data *****/

typedef struct {
  int tid;
  int start, end;
} Thread;

/***** Struct used for Nodes data *****/

typedef struct {
  double p;
  double e;
  int *From_id;
  int con_size;
  int from_size;
  int original_id;  // Store original node ID
} Node;

#define maxThreads 64

/******************** Defines ****************/
// Number of nodes
int N, num_threads;

// Convergence threashold and algorithm's parameter d
double threshold, d;

// Table of threads
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

double max_error = 1;
double sum = 0;

/***** Memory allocation - Initializations for Threads *****/

void Threads_Allocation() {

  int i;
  double N_split = (double)N / num_threads;

  // Allocate memory for threads
  Threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));

  // Stores thread's data
  Threads_data = (Thread *)malloc(num_threads * sizeof(Thread));

  // Split dataset into subsets, given to each thread
  Threads_data[0].tid = 0;
  Threads_data[0].start = 0;
  Threads_data[0].end = floor(N_split);

  for (i = 1; i < num_threads; i++) {
    Threads_data[i].tid = i;
    Threads_data[i].start = Threads_data[i - 1].end;
    if (i < (num_threads - 1)) {
      Threads_data[i].end = Threads_data[i].start + floor(N_split);
    } else {
      Threads_data[i].end = N;
    }
  }

  printf("\n");

  for (i = 0; i < num_threads; i++) {
    printf("Thread %d, start = %d, end = %d\n", Threads_data[i].tid,
           Threads_data[i].start, Threads_data[i].end);
  }

  printf("\n");
}

/***** Memory allocation - Initializations for Nodes *****/
void Nodes_Allocation() {

  int i;
  Nodes = (Node *)malloc(N * sizeof(Node));

  for (i = 0; i < N; i++) {
    Nodes[i].con_size = 0;
    Nodes[i].from_size = 0;
    Nodes[i].From_id = (int *)malloc(sizeof(int));
  }
}

/***** Read graph connections from txt file *****/
void Read_from_txt_file(char *filename) {
  FILE *fid;
  int from_idx, to_idx;
  int temp_size;
  char line[1000];
  int max_id = 0;
  int min_id = INT_MAX;

  // First pass: find min and max IDs
  fid = fopen(filename, "r");
  if (fid == NULL) {
    printf("Error opening the file\n");
    return;
  }

  while (fgets(line, sizeof(line), fid) != NULL) {
    if (sscanf(line, "%d\t%d\n", &from_idx, &to_idx)) {
      if (from_idx > max_id) max_id = from_idx;
      if (to_idx > max_id) max_id = to_idx;
      if (from_idx < min_id) min_id = from_idx;
      if (to_idx < min_id) min_id = to_idx;
    }
  }
  fclose(fid);

  // Create ID mapping
  int id_range = max_id - min_id + 1;
  int *id_to_idx = (int *)malloc(id_range * sizeof(int));
  int *idx_to_id = (int *)malloc(N * sizeof(int));
  int current_idx = 0;

  // Initialize mapping
  for (int i = 0; i < id_range; i++) {
    id_to_idx[i] = -1;  // -1 indicates unmapped
  }

  // Second pass: build mapping and read edges
  fid = fopen(filename, "r");
  if (fid == NULL) {
    printf("Error opening the file\n");
    return;
  }

  while (fgets(line, sizeof(line), fid) != NULL) {
    if (sscanf(line, "%d\t%d\n", &from_idx, &to_idx)) {
      // Map source node if not mapped
      if (id_to_idx[from_idx - min_id] == -1) {
        id_to_idx[from_idx - min_id] = current_idx;
        idx_to_id[current_idx] = from_idx;
        current_idx++;
      }
      // Map destination node if not mapped
      if (id_to_idx[to_idx - min_id] == -1) {
        id_to_idx[to_idx - min_id] = current_idx;
        idx_to_id[current_idx] = to_idx;
        current_idx++;
      }

      // Get mapped indices
      int mapped_from = id_to_idx[from_idx - min_id];
      int mapped_to = id_to_idx[to_idx - min_id];

      // Update graph structure
      Nodes[mapped_from].con_size++;
      Nodes[mapped_to].from_size++;
      temp_size = Nodes[mapped_to].from_size;
      Nodes[mapped_to].From_id = (int *)realloc(Nodes[mapped_to].From_id, temp_size * sizeof(int));
      Nodes[mapped_to].From_id[temp_size - 1] = mapped_from;
    }
  }

  printf("End of connections insertion!\n");
  fclose(fid);

  // Store original IDs
  for (int i = 0; i < N; i++) {
    Nodes[i].original_id = idx_to_id[i];
  }

  free(id_to_idx);
  free(idx_to_id);
}

/***** Read P vector from txt file*****/

void Read_P_from_txt_file() {

  FILE *fid;
  double temp_P;
  int index = 0;

  fid = fopen("P.txt", "r");
  if (fid == NULL) {
    printf("Error opening the Probabilities file\n");
  }

  while (!feof(fid)) {
    // P's values are double!
    if (fscanf(fid, "%lf\n", &temp_P)) {
      Nodes[index].p = temp_P;
      index++;
    }
  }
  printf("End of P insertion!");

  fclose(fid);
}

/***** Read E vector from txt file*****/

void Read_E_from_txt_file() {

  FILE *fid;
  double temp_E;
  int index = 0;

  fid = fopen("E.txt", "r");
  if (fid == NULL) {
    printf("Error opening the E file\n");
  }

  while (!feof(fid)) {
    // E's values are double!
    if (fscanf(fid, "%lf\n", &temp_E)) {
      Nodes[index].e = temp_E;
      index++;
    }
  }
  printf("End of E insertion!");

  fclose(fid);
}

/***** Create P and E with equal probability *****/

void Random_P_E() {

  int i;
  // Sum of P (it must be =1)
  double sum_P_1 = 0;
  // Sum of E (it must be =1)
  double sum_E_1 = 0;

  // Arrays initialization
  for (i = 0; i < N; i++) {
    Nodes[i].p = 1.0 / N;
    sum_P_1 = sum_P_1 + Nodes[i].p;

    Nodes[i].e = 1.0 / N;
    sum_E_1 = sum_E_1 + Nodes[i].e;
  }

  // Assert sum of probabilities is =1

  // Print sum of P (it must be =1)
  // printf("Sum of P = %f\n",sum_P_1);

  // Exit if sum of P is !=1
  assert(sum_P_1 = 1);

  // printf("\n");

  // Print sum of E (it must be =1)
  // printf("Sum of E = %f\n",sum_E_1);

  // Exit if sum of Pt0 is !=1
  assert(sum_E_1 = 1);
}

void *Compute_Dangling_Sum(void *arg) {
  Thread *thread_data = (Thread *)arg;
  int i;
  double local_sum = 0.0;

  for (i = thread_data->start; i < thread_data->end; i++) {
    if (Nodes[i].con_size == 0) {
      local_sum += Nodes[i].p;
    }
  }

  pthread_mutex_lock(&locksum);
  sum += local_sum;
  pthread_mutex_unlock(&locksum);

  return NULL;
}

/***** Main parallel algorithm *****/

void *Pagerank_Parallel(void *arg) {
  Thread *thread_data = (Thread *)arg;
  int i, j, from_idx;
  double local_max = 0.0;
  double rank_sum, p_old;

  for (i = thread_data->start; i < thread_data->end; i++) {
    rank_sum = 0.0;

    for (j = 0; j < Nodes[i].from_size; j++) {
      from_idx = Nodes[i].From_id[j];
      if (Nodes[from_idx].con_size > 0) {
        rank_sum += Nodes[from_idx].p / Nodes[from_idx].con_size;
      }
    }

    p_old = Nodes[i].p;

    // Handle dangling nodes: sum of ranks from them gets added outside
    Nodes[i].p = d * (rank_sum + sum) + (1.0 - d) * Nodes[i].e;

    double delta = fabs(Nodes[i].p - p_old);
    if (delta > local_max)
      local_max = delta;
  }

  // Thread-local max → global max_error
  pthread_mutex_lock(&lockmax);
  if (local_max > max_error)
    max_error = local_max;
  pthread_mutex_unlock(&lockmax);

  return NULL;
}

/***** Pagerank main algortihm *****/
void Pagerank() {
  int i;

  while (max_error > threshold) {
    iterations++;
    max_error = 0.0;
    sum = 0.0;

    // Step 1: Compute dangling node contribution
    for (i = 0; i < num_threads; i++)
      pthread_create(&Threads[i], NULL, &Compute_Dangling_Sum,
                     (void *)&Threads_data[i]);
    for (i = 0; i < num_threads; i++)
      pthread_join(Threads[i], NULL);

    sum /= N;

    // Step 2: Gauss-Seidel update
    for (i = 0; i < num_threads; i++)
      pthread_create(&Threads[i], NULL, &Pagerank_Parallel,
                     (void *)&Threads_data[i]);
    for (i = 0; i < num_threads; i++)
      pthread_join(Threads[i], NULL);

    // Step 3: Normalize PageRank vector
    double total = 0.0;
    for (i = 0; i < N; i++) {
      total += Nodes[i].p;
    }
    for (i = 0; i < N; i++) {
      Nodes[i].p /= total;
    }

    printf("Max Error in iteration %d = %f\n", iterations, max_error);
  }
}

/***** main function *****/

int main(int argc, char **argv) {

  struct timeval start, end;

  int i, j, k;
  double totaltime;

  // Check input arguments
  if (argc < 5) {
    printf("Error in arguments! Three arguments required: graph filename, N, "
           "threshold and d\n");
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
  if ((num_threads < 1) || (num_threads > maxThreads)) {
    printf("Threads number must be >= 1 and  <= %d!\n", maxThreads);
    exit(1);
  }

  Threads_Allocation();
  Nodes_Allocation();

  // OR read probabilities from files
  Read_from_txt_file(filename);
  // Read_P_from_txt_file();
  // Read_E_from_txt_file();

  Random_P_E();

  printf("\n");

  printf("Parallel version of Pagerank\n");

  gettimeofday(&start, NULL);
  Pagerank();
  gettimeofday(&end, NULL);

  // Print top 10 pages by rank
  printf("\nTop 10 pages by rank:\n");

  // Create array of indices for sorting
  int *indices = (int *)malloc(N * sizeof(int));
  for (i = 0; i < N; i++) {
    indices[i] = i;
  }

  // Sort indices based on page rank values
  for (i = 0; i < 10; i++) {
    for (int j = i + 1; j < N; j++) {
      if (Nodes[indices[j]].p > Nodes[indices[i]].p) {
        int temp = indices[i];
        indices[i] = indices[j];
        indices[j] = temp;
      }
    }
    printf("%d. Page %d (rank: %f)\n", i + 1, Nodes[indices[i]].original_id, Nodes[indices[i]].p);
  }

  free(indices);

  /*for (i = 0; i < N; i++)
  {
      printf("P_t1[%d] = %f\n",i, Nodes[i].p_t1);
  }

  printf("\n");*/

  // Print no of iterations
  printf("Total iterations: %d\n", iterations);

  totaltime =
      (((end.tv_usec - start.tv_usec) / 1.0e6 + end.tv_sec - start.tv_sec) *
       1000) /
      1000;

  printf("\nTotaltime = %f seconds\n", totaltime);

  printf("End of program!\n");

  return (EXIT_SUCCESS);
}
