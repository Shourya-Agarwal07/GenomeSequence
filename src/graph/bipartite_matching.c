#include "genome/graph/bipartite_matching.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int augment_left_vertex(size_t left,
                               const size_t *heads,
                               const size_t *neighbors,
                               const size_t *next_edges,
                               size_t *right_to_left,
                               size_t *left_to_right,
                               unsigned char *seen_right)
{
    size_t edge_index;

    for (edge_index = heads[left]; edge_index != SIZE_MAX;
         edge_index = next_edges[edge_index]) {
        size_t right = neighbors[edge_index];
        size_t previous_left;

        if (seen_right[right] != 0) {
            continue;
        }
        seen_right[right] = 1;
        previous_left = right_to_left[right];
        if (previous_left == SIZE_MAX ||
            augment_left_vertex(previous_left, heads, neighbors, next_edges,
                                right_to_left, left_to_right, seen_right)) {
            right_to_left[right] = left;
            left_to_right[left] = right;
            return 1;
        }
    }
    return 0;
}

void bipartite_matching_result_destroy(BipartiteMatchingResult *result)
{
    if (result == NULL) {
        return;
    }
    free(result->pairs);
    free(result->left_to_right);
    free(result->right_to_left);
    result->left_vertex_count = 0;
    result->right_vertex_count = 0;
    result->matching_size = 0;
    result->pairs = NULL;
    result->left_to_right = NULL;
    result->right_to_left = NULL;
}

int bipartite_maximum_matching(const BipartiteGraph *graph,
                               BipartiteMatchingResult *result)
{
    size_t *heads = NULL;
    size_t *neighbors = NULL;
    size_t *next_edges = NULL;
    size_t *left_to_right = NULL;
    size_t *right_to_left = NULL;
    unsigned char *seen_right = NULL;
    BipartiteMatch *pairs = NULL;
    size_t edge_count;
    size_t left_count;
    size_t right_count;
    size_t matching_size = 0;
    size_t i;
    int success = 0;

    if (graph == NULL || result == NULL ||
        (graph->edge_count != 0 && graph->edges == NULL)) {
        return 0;
    }
    left_count = graph->left_vertex_count;
    right_count = graph->right_vertex_count;
    edge_count = graph->edge_count;
    if (left_count == SIZE_MAX ||
        left_count > SIZE_MAX / sizeof(*heads) ||
        left_count > SIZE_MAX / sizeof(*left_to_right) ||
        right_count > SIZE_MAX / sizeof(*right_to_left) ||
        right_count > SIZE_MAX / sizeof(*seen_right) ||
        edge_count > SIZE_MAX / sizeof(*neighbors) ||
        edge_count > SIZE_MAX / sizeof(*next_edges)) {
        return 0;
    }
    for (i = 0; i < edge_count; ++i) {
        if (graph->edges[i].left >= left_count ||
            graph->edges[i].right >= right_count) {
            return 0;
        }
    }

    if (left_count != 0) {
        heads = (size_t *)malloc(left_count * sizeof(*heads));
    }
    if (left_count != 0) {
        left_to_right = (size_t *)malloc(
            left_count * sizeof(*left_to_right));
    }
    if (right_count != 0) {
        right_to_left = (size_t *)malloc(
            right_count * sizeof(*right_to_left));
        seen_right = (unsigned char *)malloc(
            right_count * sizeof(*seen_right));
    }
    if (edge_count != 0) {
        neighbors = (size_t *)malloc(edge_count * sizeof(*neighbors));
        next_edges = (size_t *)malloc(edge_count * sizeof(*next_edges));
    }
    if ((left_count != 0 && heads == NULL) ||
        (left_count != 0 && left_to_right == NULL) ||
        (right_count != 0 &&
         (right_to_left == NULL || seen_right == NULL)) ||
        (edge_count != 0 && (neighbors == NULL || next_edges == NULL))) {
        goto cleanup;
    }

    for (i = 0; i < left_count; ++i) {
        heads[i] = SIZE_MAX;
        left_to_right[i] = SIZE_MAX;
    }
    for (i = 0; i < right_count; ++i) {
        right_to_left[i] = SIZE_MAX;
    }
    for (i = 0; i < edge_count; ++i) {
        size_t left = graph->edges[i].left;

        neighbors[i] = graph->edges[i].right;
        next_edges[i] = heads[left];
        heads[left] = i;
    }

    for (i = 0; i < left_count; ++i) {
        if (right_count != 0) {
            memset(seen_right, 0, right_count * sizeof(*seen_right));
        }
        if (augment_left_vertex(i, heads, neighbors, next_edges,
                                right_to_left, left_to_right, seen_right)) {
            ++matching_size;
        }
    }

    if (matching_size != 0) {
        if (matching_size > SIZE_MAX / sizeof(*pairs)) {
            goto cleanup;
        }
        pairs = (BipartiteMatch *)malloc(
            matching_size * sizeof(*pairs));
        if (pairs == NULL) {
            goto cleanup;
        }
    }
    {
        size_t pair_index = 0;
        BipartiteMatchingResult completed;

        for (i = 0; i < left_count; ++i) {
            if (left_to_right[i] != SIZE_MAX) {
                pairs[pair_index].left = i;
                pairs[pair_index].right = left_to_right[i];
                ++pair_index;
            }
        }
        if (pair_index != matching_size) {
            goto cleanup;
        }

        completed.left_vertex_count = left_count;
        completed.right_vertex_count = right_count;
        completed.matching_size = matching_size;
        completed.pairs = pairs;
        completed.left_to_right = left_to_right;
        completed.right_to_left = right_to_left;
        pairs = NULL;
        left_to_right = NULL;
        right_to_left = NULL;
        bipartite_matching_result_destroy(result);
        *result = completed;
        success = 1;
    }

cleanup:
    free(pairs);
    free(seen_right);
    free(right_to_left);
    free(left_to_right);
    free(next_edges);
    free(neighbors);
    free(heads);
    return success;
}
