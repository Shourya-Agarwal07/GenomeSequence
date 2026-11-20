#ifndef GENOME_ALGORITHMS_DYNAMIC_PROGRAMMING_H
#define GENOME_ALGORITHMS_DYNAMIC_PROGRAMMING_H

#include <stddef.h>
#include <stdint.h>

int dp_knapsack_01(const size_t *weights,
                   const int64_t *values,
                   size_t item_count,
                   size_t capacity,
                   int64_t *maximum_value);

int dp_knapsack_01_with_selection(const size_t *weights,
                                  const int64_t *values,
                                  size_t item_count,
                                  size_t capacity,
                                  int64_t *maximum_value,
                                  size_t **selected_indices,
                                  size_t *selected_count);

int dp_longest_common_subsequence(const char *first,
                                  const char *second,
                                  size_t *lcs_length,
                                  char **lcs);

int dp_longest_common_subsequence_n(const char *first,
                                    size_t first_length,
                                    const char *second,
                                    size_t second_length,
                                    size_t *lcs_length,
                                    char **lcs);

int dp_matrix_chain_multiplication(const size_t *dimensions,
                                   size_t dimension_count,
                                   uint64_t *minimum_cost,
                                   char **parenthesization);

typedef struct DpOptimalBstNode {
    int key;
    uint64_t frequency;
    struct DpOptimalBstNode *left;
    struct DpOptimalBstNode *right;
} DpOptimalBstNode;

int dp_optimal_bst(const int *keys,
                   const uint64_t *frequencies,
                   size_t key_count,
                   uint64_t *minimum_cost,
                   DpOptimalBstNode **root);

void dp_optimal_bst_destroy(DpOptimalBstNode *root);

typedef struct {
    size_t from;
    size_t to;
    int64_t weight;
} DpBellmanFordEdge;

typedef struct {
    size_t vertex_count;
    size_t source;
    int64_t *distances;
    size_t *predecessors;
    unsigned char *reachable;
    int negative_cycle_detected;
} DpBellmanFordResult;

int dp_bellman_ford(size_t vertex_count,
                    const DpBellmanFordEdge *edges,
                    size_t edge_count,
                    size_t source,
                    DpBellmanFordResult *result);
void dp_bellman_ford_result_destroy(DpBellmanFordResult *result);

int dp_bellman_ford_reconstruct_path(const DpBellmanFordResult *result,
                                     size_t destination,
                                     size_t **path,
                                     size_t *path_length);

typedef struct {
    size_t vertex_count;
    const int64_t *weights;
    const unsigned char *edge_present;
} DpWeightedDirectedMatrix;

typedef struct {
    size_t vertex_count;
    int64_t *distances;
    size_t *next_hops;
    unsigned char *reachable;
    unsigned char *negative_cycle_affected;
    int negative_cycle_detected;
} DpFloydWarshallResult;

int dp_floyd_warshall(const DpWeightedDirectedMatrix *graph,
                      DpFloydWarshallResult *result);
void dp_floyd_warshall_result_destroy(DpFloydWarshallResult *result);

int dp_floyd_warshall_reconstruct_path(
    const DpFloydWarshallResult *result,
    size_t source,
    size_t destination,
    size_t **path,
    size_t *path_length);

#ifndef DP_TSP_MAX_VERTICES
#define DP_TSP_MAX_VERTICES 18
#endif

#if DP_TSP_MAX_VERTICES < 1 || DP_TSP_MAX_VERTICES > 20
#error DP_TSP_MAX_VERTICES must be between 1 and 20
#endif

typedef struct {
    size_t vertex_count;
    int64_t minimum_cost;
    size_t *tour;
    size_t tour_length;
} DpTspResult;

typedef enum {
    DP_TSP_ERROR = 0,
    DP_TSP_NO_TOUR = 1,
    DP_TSP_TOUR_FOUND = 2
} DpTspStatus;

DpTspStatus dp_tsp_bitmask(const DpWeightedDirectedMatrix *graph,
                           size_t start,
                           DpTspResult *result);
void dp_tsp_result_destroy(DpTspResult *result);

#endif
