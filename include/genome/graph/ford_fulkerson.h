#ifndef GENOME_GRAPH_FORD_FULKERSON_H
#define GENOME_GRAPH_FORD_FULKERSON_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t from;
    size_t to;
    int64_t capacity;
} FordFulkersonEdge;

typedef struct {
    size_t vertex_count;
    size_t edge_count;
    const FordFulkersonEdge *edges;
} FordFulkersonNetwork;

typedef struct {
    size_t from;
    size_t to;
    int64_t capacity;
} FordFulkersonCutEdge;

typedef struct {
    size_t vertex_count;
    size_t source;
    size_t sink;
    size_t edge_count;
    int64_t max_flow;
    int64_t cut_capacity;
    unsigned char *source_side;
    int64_t *edge_flows;
    FordFulkersonCutEdge *cut_edges;
    size_t cut_edge_count;
} FordFulkersonResult;

int ford_fulkerson_max_flow(const FordFulkersonNetwork *network,
                            size_t source,
                            size_t sink,
                            FordFulkersonResult *result);
void ford_fulkerson_result_destroy(FordFulkersonResult *result);

#endif
