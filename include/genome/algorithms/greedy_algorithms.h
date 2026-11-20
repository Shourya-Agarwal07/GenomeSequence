#ifndef GENOME_ALGORITHMS_GREEDY_ALGORITHMS_H
#define GENOME_ALGORITHMS_GREEDY_ALGORITHMS_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    int64_t start_time;
    int64_t finish_time;
    size_t identifier;
} Activity;

typedef struct {
    double weight;
    double value;
    size_t identifier;
} KnapsackItem;

typedef struct {
    size_t identifier;
    double amount_taken;
} KnapsackSelection;

typedef struct HuffmanNode {
    unsigned char character;
    double frequency;
    struct HuffmanNode *left;
    struct HuffmanNode *right;
} HuffmanNode;

typedef struct {
    unsigned char character;
    char *code;
} HuffmanCode;

typedef struct {
    size_t from;
    size_t to;
    double weight;
} WeightedEdge;

typedef struct {
    size_t vertex_count;
    double *adjacency;
} WeightedGraph;

typedef struct {
    size_t vertex_count;
    double *distances;
    ptrdiff_t *predecessors;
} ShortestPathResult;

typedef struct {
    size_t edge_count;
    WeightedEdge *edges;
    double total_weight;
} PrimResult;

typedef struct {
    size_t edge_count;
    WeightedEdge *edges;
    double total_weight;
    int connected;
} KruskalResult;

int greedy_activity_selection(const Activity *activities,
                              size_t count,
                              Activity *selected,
                              size_t selected_capacity,
                              size_t *selected_count);

int greedy_fractional_knapsack(const KnapsackItem *items,
                               size_t count,
                               double capacity,
                               double *selected_weights,
                               double *total_value);

HuffmanNode *greedy_huffman_build(const unsigned char *symbols,
                                  const double *frequencies,
                                  size_t count);
void greedy_huffman_destroy(HuffmanNode *root);

int greedy_huffman_encode(const unsigned char *input,
                          size_t input_length,
                          HuffmanNode *root,
                          unsigned char **encoded,
                          size_t *encoded_length);
int greedy_huffman_decode(const unsigned char *bits,
                          size_t bit_length,
                          HuffmanNode *root,
                          unsigned char **decoded,
                          size_t *decoded_length);

int greedy_dijkstra_matrix(const WeightedGraph *graph,
                           size_t source,
                           ShortestPathResult *result);
void greedy_shortest_path_destroy(ShortestPathResult *result);

int greedy_prim_matrix(const WeightedGraph *graph,
                       size_t start,
                       PrimResult *result);
void greedy_prim_destroy(PrimResult *result);

int greedy_kruskal(size_t vertex_count,
                   const WeightedEdge *edges,
                   size_t edge_count,
                   KruskalResult *result);
void greedy_kruskal_destroy(KruskalResult *result);

int activity_selection(const Activity *activities,
                       size_t count,
                       Activity *selected,
                       size_t capacity,
                       size_t *selected_count);
int fractional_knapsack(const KnapsackItem *items,
                        size_t count,
                        double capacity,
                        double *selected_weights,
                        double *total_value);
HuffmanNode *huffman_build(const unsigned char *symbols,
                           const double *frequencies,
                           size_t count);
void huffman_destroy(HuffmanNode *root);

#endif
