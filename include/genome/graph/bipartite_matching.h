#ifndef GENOME_GRAPH_BIPARTITE_MATCHING_H
#define GENOME_GRAPH_BIPARTITE_MATCHING_H

#include <stddef.h>

typedef struct {
    size_t left;
    size_t right;
} BipartiteEdge;

typedef struct {
    size_t left_vertex_count;
    size_t right_vertex_count;
    size_t edge_count;
    const BipartiteEdge *edges;
} BipartiteGraph;

typedef struct {
    size_t left;
    size_t right;
} BipartiteMatch;

typedef struct {
    size_t left_vertex_count;
    size_t right_vertex_count;
    size_t matching_size;
    BipartiteMatch *pairs;
    size_t *left_to_right;
    size_t *right_to_left;
} BipartiteMatchingResult;

int bipartite_maximum_matching(const BipartiteGraph *graph,
                               BipartiteMatchingResult *result);
void bipartite_matching_result_destroy(BipartiteMatchingResult *result);

#endif
