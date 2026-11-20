#include "genome/algorithms/greedy_algorithms.h"

#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    HuffmanNode *node;
} HuffmanHeapEntry;

enum { HUFFMAN_MAX_SYMBOLS = 256 };

static void reset_activity_selection_stats(size_t *count) {
    if (count != NULL) *count = 0;
}

static int compare_activity(const void *left_value, const void *right_value) {
    const Activity *left = (const Activity *)left_value;
    const Activity *right = (const Activity *)right_value;
    if (left->finish_time < right->finish_time) return -1;
    if (left->finish_time > right->finish_time) return 1;
    if (left->start_time < right->start_time) return -1;
    if (left->start_time > right->start_time) return 1;
    if (left->identifier < right->identifier) return -1;
    if (left->identifier > right->identifier) return 1;
    return 0;
}

int greedy_activity_selection(const Activity *activities,
                              size_t count,
                              Activity *selected,
                              size_t selected_capacity,
                              size_t *selected_count) {
    Activity *ordered;
    size_t i, used = 0;
    int64_t last_finish = 0;

    reset_activity_selection_stats(selected_count);
    if (count == 0) return 1;
    if (activities == NULL) return 0;
    if (selected == NULL && selected_capacity != 0) return 0;
    if (count > SIZE_MAX / sizeof(*ordered)) return 0;

    ordered = (Activity *)malloc(count * sizeof(*ordered));
    if (ordered == NULL) return 0;
    memcpy(ordered, activities, count * sizeof(*ordered));

    for (i = 0; i < count; ++i) {
        if (ordered[i].finish_time <= ordered[i].start_time) {
            free(ordered);
            return 0;
        }
    }
    qsort(ordered, count, sizeof(*ordered), compare_activity);

    if (selected == NULL) {
        selected_capacity = count;
    }
    for (i = 0; i < count; ++i) {
        int first = (used == 0);
        if (first || ordered[i].start_time >= last_finish) {
            if (used < selected_capacity) {
                if (selected != NULL) selected[used] = ordered[i];
            }
            ++used;
            last_finish = ordered[i].finish_time;
        }
    }

    free(ordered);
    if (selected_count != NULL) *selected_count = used;
    return used <= selected_capacity;
}

typedef struct {
    const KnapsackItem *item;
    size_t original_index;
    double ratio;
} KnapsackSortItem;

static int compare_knapsack(const void *left_value, const void *right_value) {
    const KnapsackSortItem *left = (const KnapsackSortItem *)left_value;
    const KnapsackSortItem *right = (const KnapsackSortItem *)right_value;
    if (left->ratio > right->ratio) return -1;
    if (left->ratio < right->ratio) return 1;
    if (left->item->identifier < right->item->identifier) return -1;
    if (left->item->identifier > right->item->identifier) return 1;
    return 0;
}

int greedy_fractional_knapsack(const KnapsackItem *items,
                               size_t count,
                               double capacity,
                               double *selected_weights,
                               double *total_value) {
    KnapsackSortItem *ordered;
    double remaining, total = 0.0;
    size_t i;

    if (total_value != NULL) *total_value = 0.0;
    if (!isfinite(capacity) || capacity < 0.0) return 0;
    if (count == 0) return 1;
    if (items == NULL) return 0;
    if (count > SIZE_MAX / sizeof(*ordered)) return 0;
    ordered = (KnapsackSortItem *)malloc(count * sizeof(*ordered));
    if (ordered == NULL) return 0;

    if (selected_weights != NULL) {
        for (i = 0; i < count; ++i) selected_weights[i] = 0.0;
    }
    for (i = 0; i < count; ++i) {
        if (!isfinite(items[i].weight) || items[i].weight < 0.0 ||
            !isfinite(items[i].value) || items[i].value < 0.0) {
            free(ordered);
            return 0;
        }
        ordered[i].item = &items[i];
        ordered[i].original_index = i;
        if (items[i].weight == 0.0) {
            ordered[i].ratio = items[i].value > 0.0 ? DBL_MAX : 0.0;
        } else {
            ordered[i].ratio = items[i].value / items[i].weight;
        }
    }
    qsort(ordered, count, sizeof(*ordered), compare_knapsack);
    remaining = capacity;

    for (i = 0; i < count; ++i) {
        const KnapsackSortItem *entry = &ordered[i];
        double amount = 0.0;
        if (entry->item->weight == 0.0) {
            if (entry->item->value > 0.0) {
                total += entry->item->value;
            }
            continue;
        }
        if (remaining <= 0.0) break;

        if (entry->item->weight <= remaining) {
            amount = entry->item->weight;
            remaining -= amount;
            total += entry->item->value;
        } else {
            amount = remaining;
            total += entry->item->value * (amount / entry->item->weight);
            remaining = 0.0;
        }
        if (selected_weights != NULL) {
            if (entry->item->weight == 0.0 &&
                entry->item->value > 0.0) {
                selected_weights[entry->original_index] = 0.0;
            } else {
                selected_weights[entry->original_index] = amount;
            }
        }
    }
    free(ordered);
    if (!isfinite(total)) return 0;
    if (total_value != NULL) *total_value = total;
    return 1;
}

static void huffman_heap_push(HuffmanHeapEntry *heap, size_t *count,
                              HuffmanNode *node) {
    size_t index = (*count)++;
    heap[index].node = node;
    while (index > 0) {
        size_t parent = (index - 1) / 2;
        if (heap[parent].node->frequency <= node->frequency) break;
        heap[index].node = heap[parent].node;
        index = parent;
    }
    heap[index].node = node;
}

static HuffmanNode *huffman_heap_pop(HuffmanHeapEntry *heap, size_t *count) {
    HuffmanNode *answer;
    HuffmanNode *last;
    size_t index, child;
    if (*count == 0) return NULL;
    answer = heap[0].node;
    last = heap[--(*count)].node;
    if (*count == 0) return answer;
    index = 0;
    for (;;) {
        child = index * 2 + 1;
        if (child >= *count) break;
        if (child + 1 < *count &&
            heap[child + 1].node->frequency <
                heap[child].node->frequency) ++child;
        if (heap[child].node->frequency >= last->frequency) break;
        heap[index].node = heap[child].node;
        index = child;
    }
    heap[index].node = last;
    return answer;
}

static void huffman_heap_destroy(HuffmanHeapEntry *heap, size_t count) {
    size_t i;
    for (i = 0; i < count; ++i) {
        greedy_huffman_destroy(heap[i].node);
    }
}

HuffmanNode *greedy_huffman_build(const unsigned char *symbols,
                                  const double *frequencies,
                                  size_t count) {
    HuffmanHeapEntry *heap;
    unsigned char seen[HUFFMAN_MAX_SYMBOLS] = {0};
    size_t heap_count = 0, i;
    if (count == 0) return NULL;
    if (symbols == NULL || frequencies == NULL ||
        count > HUFFMAN_MAX_SYMBOLS ||
        count > SIZE_MAX / sizeof(*heap)) return NULL;
    heap = (HuffmanHeapEntry *)malloc(count * sizeof(*heap));
    if (heap == NULL) return NULL;

    for (i = 0; i < count; ++i) {
        HuffmanNode *node;
        if (seen[symbols[i]] || !isfinite(frequencies[i]) ||
            frequencies[i] < 0.0) {
            huffman_heap_destroy(heap, heap_count);
            free(heap);
            return NULL;
        }
        seen[symbols[i]] = 1;
        node = (HuffmanNode *)calloc(1, sizeof(*node));
        if (node == NULL) {
            huffman_heap_destroy(heap, heap_count);
            free(heap);
            return NULL;
        }
        node->character = symbols[i];
        node->frequency = frequencies[i];
        huffman_heap_push(heap, &heap_count, node);
    }

    while (heap_count > 1) {
        HuffmanNode *left = huffman_heap_pop(heap, &heap_count);
        HuffmanNode *right = huffman_heap_pop(heap, &heap_count);
        HuffmanNode *parent = (HuffmanNode *)calloc(1, sizeof(*parent));
        if (parent == NULL) {
            greedy_huffman_destroy(left);
            greedy_huffman_destroy(right);
            huffman_heap_destroy(heap, heap_count);
            free(heap);
            return NULL;
        }
        parent->frequency = left->frequency + right->frequency;
        parent->left = left;
        parent->right = right;
        if (!isfinite(parent->frequency)) {
            greedy_huffman_destroy(parent);
            huffman_heap_destroy(heap, heap_count);
            free(heap);
            return NULL;
        }
        huffman_heap_push(heap, &heap_count, parent);
    }
    {
        HuffmanNode *root = huffman_heap_pop(heap, &heap_count);
        free(heap);
        return root;
    }
}

void greedy_huffman_destroy(HuffmanNode *root) {
    if (root == NULL) return;
    greedy_huffman_destroy(root->left);
    greedy_huffman_destroy(root->right);
    free(root);
}

static int huffman_encode_node(const HuffmanNode *node,
                               char *path,
                               size_t depth,
                               unsigned char **encoded,
                               size_t *length,
                               size_t *capacity,
                               unsigned char symbol) {
    size_t i;
    if (node == NULL) return 0;
    if (node->left == NULL && node->right == NULL) {
        if (node->character != symbol) return 1;
        for (i = 0; i < depth; ++i) {
            if (*length >= *capacity) {
                size_t new_capacity = *capacity == 0 ? 16 : *capacity * 2;
                unsigned char *larger;
                if (new_capacity < *capacity ||
                    new_capacity > SIZE_MAX / sizeof(**encoded)) return 0;
                larger = (unsigned char *)realloc(*encoded,
                                                  new_capacity * sizeof(**encoded));
                if (larger == NULL) return 0;
                *encoded = larger;
                *capacity = new_capacity;
            }
            (*encoded)[(*length)++] = (unsigned char)(path[i] == '1');
        }
        if (depth == 0) {
            if (*length >= *capacity) {
                unsigned char *larger =
                    (unsigned char *)realloc(*encoded, 1);
                if (larger == NULL) return 0;
                *encoded = larger;
                *capacity = 1;
            }
            (*encoded)[(*length)++] = 0;
        }
        return 1;
    }
    if (depth >= HUFFMAN_MAX_SYMBOLS - 1) return 0;
    path[depth] = '0';
    if (!huffman_encode_node(node->left, path, depth + 1, encoded,
                             length, capacity, symbol)) return 0;
    path[depth] = '1';
    return huffman_encode_node(node->right, path, depth + 1, encoded,
                               length, capacity, symbol);
}

int greedy_huffman_encode(const unsigned char *input,
                          size_t input_length,
                          HuffmanNode *root,
                          unsigned char **encoded,
                          size_t *encoded_length) {
    size_t i;
    char path[256];
    unsigned char *all = NULL;
    size_t length = 0, capacity = 0;
    if (encoded != NULL) *encoded = NULL;
    if (encoded_length != NULL) *encoded_length = 0;
    if (root == NULL) return 0;
    if (input_length != 0 && input == NULL) return 0;
    if (encoded == NULL || encoded_length == NULL) return 0;

    for (i = 0; i < input_length; ++i) {
        size_t start = length;
        if (!huffman_encode_node(root, path, 0, &all, &length, &capacity,
                                 input[i])) {
            free(all);
            return 0;
        }
        if (length == start &&
            !(root->left == NULL && root->right == NULL &&
              root->character == input[i])) {
            free(all);
            return 0;
        }
    }
    if (length == 0 && input_length == 0) {
        all = (unsigned char *)malloc(1);
        if (all == NULL) return 0;
        all[0] = 0;
    }
    *encoded = all;
    *encoded_length = length;
    return 1;
}

int greedy_huffman_decode(const unsigned char *bits,
                          size_t bit_length,
                          HuffmanNode *root,
                          unsigned char **decoded,
                          size_t *decoded_length) {
    HuffmanNode *node;
    unsigned char *output = NULL;
    size_t length = 0, capacity = 0, i;
    if (decoded != NULL) *decoded = NULL;
    if (decoded_length != NULL) *decoded_length = 0;
    if (root == NULL || decoded == NULL || decoded_length == NULL) return 0;
    if (bit_length != 0 && bits == NULL) return 0;
    node = root;

    if (root->left == NULL && root->right == NULL) {
        if (root->left != NULL || root->right != NULL) return 0;
        if (bit_length == SIZE_MAX) return 0;
        output = (unsigned char *)malloc(bit_length + 1);
        if (output == NULL) return 0;
        for (i = 0; i < bit_length; ++i) {
            if (bits[i] != 0) {
                free(output);
                return 0;
            }
            output[i] = root->character;
        }
        output[bit_length] = 0;
        *decoded = output;
        *decoded_length = bit_length;
        return 1;
    }
    for (i = 0; i < bit_length; ++i) {
        if (bits[i] > 1) {
            free(output);
            return 0;
        }
        node = bits[i] ? node->right : node->left;
        if (node == NULL) {
            free(output);
            return 0;
        }
        if (node->left == NULL && node->right == NULL) {
            if (length + 1 >= capacity) {
                size_t next = capacity == 0 ? 16 : capacity * 2;
                unsigned char *larger;
                if (next < capacity || next > SIZE_MAX - 1) {
                    free(output); return 0;
                }
                larger = (unsigned char *)realloc(output, next);
                if (larger == NULL) { free(output); return 0; }
                output = larger;
                capacity = next;
            }
            output[length++] = node->character;
            node = root;
        }
    }
    if (node != root) { free(output); return 0; }
    if (length + 1 < length) { free(output); return 0; }
    if (length == 0) {
        unsigned char *larger = (unsigned char *)realloc(output, 1);
        if (larger == NULL) {
            free(output);
            return 0;
        }
        output = larger;
    } else {
        unsigned char *larger = (unsigned char *)realloc(output, length + 1);
        if (larger == NULL) { free(output); return 0; }
        output = larger;
    }
    output[length] = 0;
    *decoded = output;
    *decoded_length = length;
    return 1;
}

static int graph_has_valid_storage(const WeightedGraph *graph) {
    size_t count;
    if (graph == NULL || graph->vertex_count == 0 ||
        graph->adjacency == NULL) {
        return 0;
    }
    count = graph->vertex_count;
    return count <= SIZE_MAX / count &&
           count <= (SIZE_MAX / sizeof(double)) / count &&
           count <= SIZE_MAX / sizeof(double) &&
           count <= SIZE_MAX / sizeof(ptrdiff_t) &&
           count - 1 <= SIZE_MAX / sizeof(WeightedEdge) &&
           (uintmax_t)(count - 1) <= (uintmax_t)PTRDIFF_MAX;
}

static int graph_is_valid_for_dijkstra(const WeightedGraph *graph) {
    size_t i, j, count;
    if (!graph_has_valid_storage(graph)) return 0;
    count = graph->vertex_count;
    for (i = 0; i < count; ++i) {
        for (j = 0; j < count; ++j) {
            double weight = graph->adjacency[i * count + j];
            if (!isfinite(weight) || weight < 0.0) return 0;
        }
    }
    return 1;
}

static int graph_is_valid_for_prim(const WeightedGraph *graph) {
    size_t i, j, count;
    if (!graph_has_valid_storage(graph)) return 0;
    count = graph->vertex_count;
    for (i = 0; i < count; ++i) {
        for (j = 0; j < count; ++j) {
            double weight = graph->adjacency[i * count + j];
            if (!isfinite(weight) ||
                weight != graph->adjacency[j * count + i]) return 0;
        }
    }
    return 1;
}

int greedy_dijkstra_matrix(const WeightedGraph *graph,
                           size_t source,
                           ShortestPathResult *result) {
    unsigned char *visited;
    size_t i, count;
    double *distances;
    ptrdiff_t *predecessors;
    ShortestPathResult completed;
    if (result == NULL) return 0;
    if (!graph_is_valid_for_dijkstra(graph) ||
        source >= graph->vertex_count) return 0;
    count = graph->vertex_count;
    distances = (double *)malloc(count * sizeof(*distances));
    predecessors = (ptrdiff_t *)malloc(count * sizeof(*predecessors));
    visited = (unsigned char *)calloc(count, sizeof(*visited));
    if (distances == NULL || predecessors == NULL || visited == NULL) {
        free(distances); free(predecessors); free(visited); return 0;
    }
    for (i = 0; i < count; ++i) {
        distances[i] = DBL_MAX;
        predecessors[i] = -1;
    }
    distances[source] = 0.0;
    for (i = 0; i < count; ++i) {
        size_t u = SIZE_MAX, j;
        for (j = 0; j < count; ++j) {
            if (!visited[j] &&
                (u == SIZE_MAX || distances[j] < distances[u])) u = j;
        }
        if (u == SIZE_MAX || distances[u] == DBL_MAX) break;
        visited[u] = 1;
        for (j = 0; j < count; ++j) {
            double weight = graph->adjacency[u * count + j];
            if (!visited[j] && isfinite(weight) && weight >= 0.0 &&
                distances[u] <= DBL_MAX - weight &&
                distances[u] + weight < distances[j]) {
                distances[j] = distances[u] + weight;
                predecessors[j] = (ptrdiff_t)u;
            }
        }
    }
    free(visited);
    completed.vertex_count = count;
    completed.distances = distances;
    completed.predecessors = predecessors;
    greedy_shortest_path_destroy(result);
    *result = completed;
    return 1;
}

void greedy_shortest_path_destroy(ShortestPathResult *result) {
    if (result == NULL) return;
    free(result->distances);
    free(result->predecessors);
    result->vertex_count = 0;
    result->distances = NULL;
    result->predecessors = NULL;
}

int greedy_prim_matrix(const WeightedGraph *graph,
                       size_t start,
                       PrimResult *result) {
    unsigned char *used;
    double *key;
    ptrdiff_t *parent;
    size_t i, count;
    size_t edge_count = 0;
    double total_weight = 0.0;
    WeightedEdge *edges = NULL;
    PrimResult completed;
    if (result == NULL) return 0;
    if (!graph_is_valid_for_prim(graph) || start >= graph->vertex_count) {
        return 0;
    }
    count = graph->vertex_count;
    used = (unsigned char *)calloc(count, sizeof(*used));
    key = (double *)malloc(count * sizeof(*key));
    parent = (ptrdiff_t *)malloc(count * sizeof(*parent));
    if (used == NULL || key == NULL || parent == NULL) {
        free(used); free(key); free(parent); return 0;
    }
    for (i = 0; i < count; ++i) {
        key[i] = DBL_MAX;
        parent[i] = -1;
    }
    key[start] = 0.0;
    for (i = 0; i < count; ++i) {
        size_t u = SIZE_MAX, j;
        for (j = 0; j < count; ++j) {
            if (!used[j] && (u == SIZE_MAX || key[j] < key[u])) u = j;
        }
        if (u == SIZE_MAX) break;
        if (key[u] == DBL_MAX) {
            key[u] = 0.0;
            parent[u] = -1;
        }
        used[u] = 1;
        for (j = 0; j < count; ++j) {
            double weight = graph->adjacency[u * count + j];
            if (!used[j] && isfinite(weight) && weight < key[j]) {
                key[j] = weight;
                parent[j] = (ptrdiff_t)u;
            }
        }
    }
    edges = (WeightedEdge *)malloc((count - 1) * sizeof(*edges));
    if (count > 1 && edges == NULL) {
        free(used); free(key); free(parent); return 0;
    }
    for (i = 0; i < count; ++i) {
        if (i != start && parent[i] >= 0) {
            edges[edge_count].from = (size_t)parent[i];
            edges[edge_count].to = i;
            edges[edge_count].weight = key[i];
            if (!isfinite(total_weight + key[i])) {
                free(edges);
                free(used); free(key); free(parent);
                return 0;
            }
            total_weight += key[i];
            ++edge_count;
        }
    }
    free(used); free(key); free(parent);
    completed.edge_count = edge_count;
    completed.edges = edges;
    completed.total_weight = total_weight;
    greedy_prim_destroy(result);
    *result = completed;
    return 1;
}

void greedy_prim_destroy(PrimResult *result) {
    if (result == NULL) return;
    free(result->edges);
    result->edge_count = 0;
    result->edges = NULL;
    result->total_weight = 0.0;
}

typedef struct {
    size_t parent;
    size_t rank;
} DisjointSetNode;

static int compare_weighted_edges(const void *left_value,
                                  const void *right_value) {
    const WeightedEdge *left = (const WeightedEdge *)left_value;
    const WeightedEdge *right = (const WeightedEdge *)right_value;
    if (left->weight < right->weight) return -1;
    if (left->weight > right->weight) return 1;
    if (left->from < right->from) return -1;
    if (left->from > right->from) return 1;
    if (left->to < right->to) return -1;
    if (left->to > right->to) return 1;
    return 0;
}

static size_t disjoint_set_find(DisjointSetNode *sets, size_t item) {
    size_t root = item;
    while (sets[root].parent != root) {
        root = sets[root].parent;
    }
    while (sets[item].parent != item) {
        size_t next = sets[item].parent;
        sets[item].parent = root;
        item = next;
    }
    return root;
}

static int disjoint_set_union(DisjointSetNode *sets,
                              size_t left,
                              size_t right) {
    size_t left_root = disjoint_set_find(sets, left);
    size_t right_root = disjoint_set_find(sets, right);
    if (left_root == right_root) return 0;
    if (sets[left_root].rank < sets[right_root].rank) {
        sets[left_root].parent = right_root;
    } else if (sets[left_root].rank > sets[right_root].rank) {
        sets[right_root].parent = left_root;
    } else {
        sets[right_root].parent = left_root;
        ++sets[left_root].rank;
    }
    return 1;
}

int greedy_kruskal(size_t vertex_count,
                   const WeightedEdge *edges,
                   size_t edge_count,
                   KruskalResult *result) {
    WeightedEdge *ordered = NULL;
    WeightedEdge *selected = NULL;
    DisjointSetNode *sets = NULL;
    size_t i, selected_count = 0, selected_capacity;
    double total_weight = 0.0;
    KruskalResult completed;

    if (result == NULL) return 0;

    if ((edge_count != 0 && edges == NULL) ||
        (vertex_count == 0 && edge_count != 0) ||
        vertex_count > SIZE_MAX / sizeof(*sets) ||
        edge_count > SIZE_MAX / sizeof(*ordered)) {
        return 0;
    }
    for (i = 0; i < edge_count; ++i) {
        if (edges[i].from >= vertex_count || edges[i].to >= vertex_count ||
            !isfinite(edges[i].weight)) {
            return 0;
        }
    }

    if (edge_count != 0) {
        ordered = (WeightedEdge *)malloc(edge_count * sizeof(*ordered));
        if (ordered == NULL) return 0;
        memcpy(ordered, edges, edge_count * sizeof(*ordered));
        if (edge_count > 1) {
            qsort(ordered, edge_count, sizeof(*ordered),
                  compare_weighted_edges);
        }
    }
    if (vertex_count != 0) {
        sets = (DisjointSetNode *)malloc(vertex_count * sizeof(*sets));
        if (sets == NULL) {
            free(ordered);
            return 0;
        }
        for (i = 0; i < vertex_count; ++i) {
            sets[i].parent = i;
            sets[i].rank = 0;
        }
    }

    selected_capacity = vertex_count == 0 ? 0 : vertex_count - 1;
    if (selected_capacity != 0) {
        if (selected_capacity > SIZE_MAX / sizeof(*selected)) {
            free(sets);
            free(ordered);
            return 0;
        }
        selected = (WeightedEdge *)malloc(selected_capacity *
                                           sizeof(*selected));
        if (selected == NULL) {
            free(sets);
            free(ordered);
            return 0;
        }
    }

    for (i = 0; i < edge_count && selected_count < selected_capacity; ++i) {
        double new_total;
        if (disjoint_set_union(sets, ordered[i].from, ordered[i].to)) {
            new_total = total_weight + ordered[i].weight;
            if (!isfinite(new_total)) {
                free(selected);
                free(sets);
                free(ordered);
                return 0;
            }
            selected[selected_count++] = ordered[i];
            total_weight = new_total;
        }
    }

    free(sets);
    free(ordered);
    completed.edge_count = selected_count;
    completed.edges = selected;
    completed.total_weight = total_weight;
    completed.connected = vertex_count != 0 &&
                          selected_count == vertex_count - 1;
    greedy_kruskal_destroy(result);
    *result = completed;
    return 1;
}

void greedy_kruskal_destroy(KruskalResult *result) {
    if (result == NULL) return;
    free(result->edges);
    result->edge_count = 0;
    result->edges = NULL;
    result->total_weight = 0.0;
    result->connected = 0;
}

int activity_selection(const Activity *activities,
                       size_t count,
                       Activity *selected,
                       size_t capacity,
                       size_t *selected_count) {
    return greedy_activity_selection(activities, count, selected,
                                     capacity, selected_count);
}

int fractional_knapsack(const KnapsackItem *items,
                        size_t count,
                        double capacity,
                        double *selected_weights,
                        double *total_value) {
    return greedy_fractional_knapsack(items, count, capacity,
                                      selected_weights, total_value);
}

HuffmanNode *huffman_build(const unsigned char *symbols,
                           const double *frequencies,
                           size_t count) {
    return greedy_huffman_build(symbols, frequencies, count);
}

void huffman_destroy(HuffmanNode *root) {
    greedy_huffman_destroy(root);
}
