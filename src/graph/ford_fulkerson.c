#include "genome/graph/ford_fulkerson.h"

#include <limits.h>
#include <stdlib.h>

typedef struct {
    size_t from;
    size_t to;
    size_t next;
    int64_t residual_capacity;
} ResidualArc;

static int augmenting_path_bfs(const ResidualArc *arcs,
                               const size_t *heads,
                               size_t vertex_count,
                               size_t source,
                               size_t sink,
                               size_t *parent_arc,
                               size_t *queue)
{
    size_t front = 0;
    size_t back = 0;
    size_t vertex;

    for (vertex = 0; vertex < vertex_count; ++vertex) {
        parent_arc[vertex] = SIZE_MAX;
    }
    parent_arc[source] = SIZE_MAX - 1;
    queue[back++] = source;

    while (front < back) {
        size_t current = queue[front++];
        size_t arc_index;

        for (arc_index = heads[current]; arc_index != SIZE_MAX;
             arc_index = arcs[arc_index].next) {
            size_t next = arcs[arc_index].to;

            if (arcs[arc_index].residual_capacity <= 0 ||
                parent_arc[next] != SIZE_MAX) {
                continue;
            }
            parent_arc[next] = arc_index;
            if (next == sink) {
                return 1;
            }
            queue[back++] = next;
        }
    }
    return 0;
}

static void ford_fulkerson_mark_source_side(const ResidualArc *arcs,
                                            const size_t *heads,
                                            size_t source,
                                            unsigned char *source_side,
                                            size_t *queue)
{
    size_t front = 0;
    size_t back = 0;

    source_side[source] = 1;
    queue[back++] = source;
    while (front < back) {
        size_t current = queue[front++];
        size_t arc_index;

        for (arc_index = heads[current]; arc_index != SIZE_MAX;
             arc_index = arcs[arc_index].next) {
            size_t next = arcs[arc_index].to;

            if (arcs[arc_index].residual_capacity > 0 &&
                source_side[next] == 0) {
                source_side[next] = 1;
                queue[back++] = next;
            }
        }
    }
}

void ford_fulkerson_result_destroy(FordFulkersonResult *result)
{
    if (result == NULL) {
        return;
    }
    free(result->source_side);
    free(result->edge_flows);
    free(result->cut_edges);
    result->vertex_count = 0;
    result->source = 0;
    result->sink = 0;
    result->edge_count = 0;
    result->max_flow = 0;
    result->cut_capacity = 0;
    result->source_side = NULL;
    result->edge_flows = NULL;
    result->cut_edges = NULL;
    result->cut_edge_count = 0;
}

int ford_fulkerson_max_flow(const FordFulkersonNetwork *network,
                            size_t source,
                            size_t sink,
                            FordFulkersonResult *result)
{
    ResidualArc *arcs = NULL;
    size_t *heads = NULL;
    size_t *parent_arc = NULL;
    size_t *queue = NULL;
    unsigned char *source_side = NULL;
    int64_t *edge_flows = NULL;
    FordFulkersonCutEdge *cut_edges = NULL;
    size_t vertex_count;
    size_t edge_count;
    size_t arc_count;
    size_t i;
    size_t cut_edge_count = 0;
    int64_t max_flow = 0;
    int64_t cut_capacity = 0;
    int success = 0;

    if (network == NULL || result == NULL ||
        network->vertex_count == 0 ||
        source >= network->vertex_count ||
        sink >= network->vertex_count || source == sink ||
        (network->edge_count != 0 && network->edges == NULL)) {
        return 0;
    }
    vertex_count = network->vertex_count;
    edge_count = network->edge_count;
    if (vertex_count > SIZE_MAX / sizeof(*heads) ||
        vertex_count > SIZE_MAX / sizeof(*parent_arc) ||
        vertex_count > SIZE_MAX / sizeof(*queue) ||
        vertex_count > SIZE_MAX / sizeof(*source_side) ||
        edge_count > SIZE_MAX / 2 ||
        edge_count > SIZE_MAX / sizeof(*edge_flows) ||
        edge_count > SIZE_MAX / sizeof(*cut_edges)) {
        return 0;
    }
    arc_count = edge_count * 2;
    if (arc_count > SIZE_MAX / sizeof(*arcs)) {
        return 0;
    }
    for (i = 0; i < edge_count; ++i) {
        if (network->edges[i].from >= vertex_count ||
            network->edges[i].to >= vertex_count ||
            network->edges[i].capacity < 0) {
            return 0;
        }
    }

    heads = (size_t *)malloc(vertex_count * sizeof(*heads));
    parent_arc = (size_t *)malloc(vertex_count * sizeof(*parent_arc));
    queue = (size_t *)malloc(vertex_count * sizeof(*queue));
    source_side = (unsigned char *)calloc(
        vertex_count, sizeof(*source_side));
    if (edge_count != 0) {
        arcs = (ResidualArc *)malloc(arc_count * sizeof(*arcs));
    }
    if (heads == NULL || parent_arc == NULL || queue == NULL ||
        source_side == NULL || (edge_count != 0 && arcs == NULL)) {
        goto cleanup;
    }

    for (i = 0; i < vertex_count; ++i) {
        heads[i] = SIZE_MAX;
    }
    for (i = 0; i < edge_count; ++i) {
        size_t forward = i * 2;
        size_t reverse = forward + 1;
        const FordFulkersonEdge *input = &network->edges[i];

        arcs[forward].from = input->from;
        arcs[forward].to = input->to;
        arcs[forward].residual_capacity = input->capacity;
        arcs[forward].next = heads[input->from];
        heads[input->from] = forward;

        arcs[reverse].from = input->to;
        arcs[reverse].to = input->from;
        arcs[reverse].residual_capacity = 0;
        arcs[reverse].next = heads[input->to];
        heads[input->to] = reverse;
    }

    while (augmenting_path_bfs(arcs, heads, vertex_count, source, sink,
                               parent_arc, queue)) {
        int64_t bottleneck = INT64_MAX;
        size_t current = sink;

        while (current != source) {
            size_t arc_index = parent_arc[current];

            if (arc_index == SIZE_MAX || arc_index == SIZE_MAX - 1) {
                goto cleanup;
            }
            if (arcs[arc_index].residual_capacity < bottleneck) {
                bottleneck = arcs[arc_index].residual_capacity;
            }
            current = arcs[arc_index].from;
        }
        if (bottleneck <= 0 || bottleneck > INT64_MAX - max_flow) {
            goto cleanup;
        }

        current = sink;
        while (current != source) {
            size_t arc_index = parent_arc[current];
            size_t reverse_index = arc_index ^ (size_t)1;

            if (bottleneck >
                INT64_MAX - arcs[reverse_index].residual_capacity) {
                goto cleanup;
            }
            arcs[arc_index].residual_capacity -= bottleneck;
            arcs[reverse_index].residual_capacity += bottleneck;
            current = arcs[arc_index].from;
        }
        max_flow += bottleneck;
    }

    ford_fulkerson_mark_source_side(arcs, heads, source,
                                    source_side, queue);
    if (edge_count != 0) {
        edge_flows = (int64_t *)malloc(edge_count * sizeof(*edge_flows));
        if (edge_flows == NULL) {
            goto cleanup;
        }
    }
    for (i = 0; i < edge_count; ++i) {
        const FordFulkersonEdge *input = &network->edges[i];

        edge_flows[i] = arcs[i * 2 + 1].residual_capacity;
        if (source_side[input->from] != 0 &&
            source_side[input->to] == 0) {
            ++cut_edge_count;
            if (input->capacity > INT64_MAX - cut_capacity) {
                goto cleanup;
            }
            cut_capacity += input->capacity;
        }
    }
    if (cut_capacity != max_flow) {
        goto cleanup;
    }
    if (cut_edge_count != 0) {
        cut_edges = (FordFulkersonCutEdge *)malloc(
            cut_edge_count * sizeof(*cut_edges));
        if (cut_edges == NULL) {
            goto cleanup;
        }
        cut_edge_count = 0;
        for (i = 0; i < edge_count; ++i) {
            const FordFulkersonEdge *input = &network->edges[i];

            if (source_side[input->from] != 0 &&
                source_side[input->to] == 0) {
                cut_edges[cut_edge_count].from = input->from;
                cut_edges[cut_edge_count].to = input->to;
                cut_edges[cut_edge_count].capacity = input->capacity;
                ++cut_edge_count;
            }
        }
    }

    {
        FordFulkersonResult completed;

        completed.vertex_count = vertex_count;
        completed.source = source;
        completed.sink = sink;
        completed.edge_count = edge_count;
        completed.max_flow = max_flow;
        completed.cut_capacity = cut_capacity;
        completed.source_side = source_side;
        completed.edge_flows = edge_flows;
        completed.cut_edges = cut_edges;
        completed.cut_edge_count = cut_edge_count;
        source_side = NULL;
        edge_flows = NULL;
        cut_edges = NULL;
        ford_fulkerson_result_destroy(result);
        *result = completed;
        success = 1;
    }

cleanup:
    free(cut_edges);
    free(edge_flows);
    free(source_side);
    free(queue);
    free(parent_arc);
    free(heads);
    free(arcs);
    return success;
}
