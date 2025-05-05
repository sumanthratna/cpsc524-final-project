// g++ -O3 pagerank_serial.cpp -o pagerank
// ./pagerank graph.txt 916428 0.000001 0.85

/********************************************************************/
/*    Pagerank project 2014 - Serial version                        */
/*    	*based on Cleve Moler's matlab implementation               */
/*                                                                  */
/*    Implemented by Nikos Katirtzis (nikos912000)                  */
/********************************************************************/

/******************** Includes - Defines ****************/
#include <assert.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

/***** Struct for timestamps *****/
struct timeval start, end;

/***** Struct used for Nodes data *****/

typedef struct {
  double p_t0;
  double p_t1;
  double e;
  int *To_id;
  int con_size;
  int original_id; // Store original node ID
} Node;

/******************** Defines ****************/
// Number of nodes
int N;

// Convergence threashold and algorithm's parameter d
double threshold, d;

// Table of node's data
Node *Nodes;

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
    printf("Error opening data file\n");
    return;
  }

  while (fgets(line, sizeof(line), fid) != NULL) {
    if (strncmp(line, "#", 1) != 0) {
      if (sscanf(line, "%d\t%d\n", &from_idx, &to_idx)) {
        if (from_idx > max_id)
          max_id = from_idx;
        if (to_idx > max_id)
          max_id = to_idx;
        if (from_idx < min_id)
          min_id = from_idx;
        if (to_idx < min_id)
          min_id = to_idx;
      }
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
    id_to_idx[i] = -1; // -1 indicates unmapped
  }

  // Second pass: build mapping and read edges
  fid = fopen(filename, "r");
  if (fid == NULL) {
    printf("Error opening data file\n");
    return;
  }

  while (fgets(line, sizeof(line), fid) != NULL) {
    if (strncmp(line, "#", 1) != 0) {
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
        temp_size = Nodes[mapped_from].con_size;
        Nodes[mapped_from].To_id =
            (int *)realloc(Nodes[mapped_from].To_id, temp_size * sizeof(int));
        Nodes[mapped_from].To_id[temp_size - 1] = mapped_to;
      }
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
      Nodes[index].p_t1 = temp_P;
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
    Nodes[i].p_t0 = 0;
    Nodes[i].p_t1 = 1;
    Nodes[i].p_t1 = (double)Nodes[i].p_t1 / N;

    sum_P_1 = sum_P_1 + Nodes[i].p_t1;

    Nodes[i].e = 1;
    Nodes[i].e = (double)Nodes[i].e / N;
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

/***** Main function *****/

int main(int argc, char **argv) {
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

  int i, j, k;
  double totaltime;

  // a constant value contributed of all nodes with connectivity = 0
  // it's going to be addes to all node's new probability
  double sum = 0;

  // Allocate memory for N nodes
  Nodes = (Node *)malloc(N * sizeof(Node));

  for (i = 0; i < N; i++) {
    Nodes[i].con_size = 0;
    Nodes[i].To_id = (int *)malloc(sizeof(int));
  }

  Read_from_txt_file(filename);

  // set random probabilities
  Random_P_E();
  // OR read probabilities from files
  // Read_P_from_txt_file();
  // Read_E_from_txt_file();

  gettimeofday(&start, NULL);

  /********** Start of algorithm **********/

  // Iterations counter
  int iterations = 0;
  int index;

  // Or any value > threshold
  double max_error = 1;

  printf("\nSerial version of Pagerank\n");

  // Continue if we don't have convergence yet
  while (max_error > threshold) {
    sum = 0;

    // Initialize P(t) and P(t + 1) values
    for (i = 0; i < N; i++) {
      // Update the "old" P table with the new one
      Nodes[i].p_t0 = Nodes[i].p_t1;
      Nodes[i].p_t1 = 0;
    }

    // Find P for each webpage
    for (i = 0; i < N; i++) {

      if (Nodes[i].con_size != 0) {

        // Compute the total probability, contributed by node's neighbors
        for (j = 0; j < Nodes[i].con_size; j++) {
          index = Nodes[i].To_id[j];
          Nodes[index].p_t1 =
              Nodes[index].p_t1 + (double)Nodes[i].p_t0 / Nodes[i].con_size;
        }

      }

      else {
        // Contribute to all
        sum = sum + (double)Nodes[i].p_t0 / N;
      }
    }

    max_error = -1;

    // Compute the new probabilities and find maximum error
    for (i = 0; i < N; i++) {
      Nodes[i].p_t1 = d * (Nodes[i].p_t1 + sum) + (1 - d) * Nodes[i].e;

      if (fabs(Nodes[i].p_t1 - Nodes[i].p_t0) > max_error) {
        max_error = fabs(Nodes[i].p_t1 - Nodes[i].p_t0);
      }
    }

    printf("Max Error in iteration %d = %f\n", iterations + 1, max_error);
    iterations++;
  }

  gettimeofday(&end, NULL);

  printf("\n");

  // Print final probabilitities
  /*for (i = 0; i < N; i++)
  {
      printf("P_t1[%d] = %f\n",i,Nodes[i].p_t1);
  }

  printf("\n");*/

  // Print no of iterations
  printf("Total iterations: %d\n", iterations);

  totaltime =
      (((end.tv_usec - start.tv_usec) / 1.0e6 + end.tv_sec - start.tv_sec) *
       1000) /
      1000;

  printf("\nTotaltime = %f seconds\n", totaltime);

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
      if (Nodes[indices[j]].p_t1 > Nodes[indices[i]].p_t1) {
        int temp = indices[i];
        indices[i] = indices[j];
        indices[j] = temp;
      }
    }
    printf("%d. Page %d (rank: %f)\n", i + 1, Nodes[indices[i]].original_id,
           Nodes[indices[i]].p_t1);
  }

  printf("End of program!\n");

  return (EXIT_SUCCESS);
}