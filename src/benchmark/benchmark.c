#include "genome/benchmark/benchmark.h"

#include "genome/algorithms/backtracking_algorithms.h"
#include "genome/algorithms/dynamic_programming.h"
#include "genome/algorithms/greedy_algorithms.h"
#include "genome/algorithms/sequence_algorithms.h"
#include "genome/graph/ford_fulkerson.h"
#include "genome/string_matching/string_matching.h"

#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SEQUENCE_BENCHMARK_IDENTIFIER_CAPACITY \
    (3U * sizeof(size_t) + sizeof("seq"))

int benchmark_measure(const char *algorithm,
                      size_t input_size,
                      size_t repetitions,
                      BenchmarkCallback prepare,
                      BenchmarkCallback operation,
                      BenchmarkCallback verify,
                      BenchmarkCleanup cleanup,
                      BenchmarkOperationCount operation_count,
                      void *context,
                      BenchmarkRecord *record)
{
    BenchmarkRecord completed = {0};
    long double total_seconds = 0.0L;
    uint64_t counted_operations = 0;
    size_t i;
    size_t name_length;

    if (algorithm == NULL || repetitions == 0 || prepare == NULL ||
        operation == NULL || verify == NULL || cleanup == NULL ||
        context == NULL || record == NULL) {
        return 0;
    }
    name_length = strlen(algorithm);
    if (name_length == 0 ||
        name_length >= sizeof(completed.algorithm)) {
        return 0;
    }

    for (i = 0; i < repetitions; ++i) {
        struct timespec start_time;
        struct timespec end_time;
        long double elapsed;
        uint64_t current_operations = 0;
        int prepared = prepare(context);
        int ran = 0;
        int ended = 0;

        if (prepared && timespec_get(&start_time, TIME_UTC) == TIME_UTC) {
            ran = operation(context);
            ended = timespec_get(&end_time, TIME_UTC) == TIME_UTC;
        }
        if (ran && ended) {
            elapsed = (long double)difftime(end_time.tv_sec,
                                             start_time.tv_sec) +
                      ((long double)end_time.tv_nsec -
                       (long double)start_time.tv_nsec) /
                          1000000000.0L;
        } else {
            elapsed = -1.0L;
        }
        if (!prepared || !ran || !ended || elapsed < 0.0L) {
            cleanup(context);
            return 0;
        }
        if (!verify(context)) {
            cleanup(context);
            return 0;
        }
        if (operation_count != NULL &&
            !operation_count(context, &current_operations)) {
            cleanup(context);
            return 0;
        }
        cleanup(context);

        total_seconds += elapsed;
        if (!isfinite(total_seconds) || total_seconds > (long double)DBL_MAX) {
            return 0;
        }
        if (operation_count != NULL) {
            if (i != 0 && current_operations != counted_operations) {
                return 0;
            }
            counted_operations = current_operations;
        }
    }

    memcpy(completed.algorithm, algorithm, name_length + 1);
    completed.input_size = input_size;
    completed.average_seconds =
        (double)(total_seconds / (long double)repetitions);
    completed.repetitions = repetitions;
    completed.verified = 1;
    completed.operation_count = counted_operations;
    completed.has_operation_count = operation_count != NULL;
    *record = completed;
    return 1;
}

typedef enum {
    SEQUENCE_BENCHMARK_MERGE,
    SEQUENCE_BENCHMARK_QUICK,
    SEQUENCE_BENCHMARK_ITERATIVE_SEARCH,
    SEQUENCE_BENCHMARK_RECURSIVE_SEARCH
} SequenceBenchmarkKind;

typedef struct {
    const DnaSequence *source;
    DnaSequence *working;
    size_t count;
    SequenceSortCriterion criterion;
    SequenceBenchmarkKind kind;
    DnaSequence target;
    SequenceSearchResult search_result;
    size_t search_index;
    SequenceAlgorithmStats statistics;
} SequenceBenchmarkContext;

static int sequence_benchmark_prepare(void *context_pointer)
{
    SequenceBenchmarkContext *context =
        (SequenceBenchmarkContext *)context_pointer;

    free(context->working);
    context->working = NULL;
    memset(&context->statistics, 0, sizeof(context->statistics));
    context->search_index = SIZE_MAX;
    context->search_result = SEQUENCE_SEARCH_ERROR;
    if (context->kind == SEQUENCE_BENCHMARK_MERGE ||
        context->kind == SEQUENCE_BENCHMARK_QUICK) {
        if (context->count > SIZE_MAX / sizeof(*context->working)) {
            return 0;
        }
        context->working = (DnaSequence *)malloc(
            context->count * sizeof(*context->working));
        if (context->working == NULL) {
            return 0;
        }
        memcpy(context->working, context->source,
               context->count * sizeof(*context->working));
    }
    return 1;
}

static int sequence_benchmark_operation(void *context_pointer)
{
    SequenceBenchmarkContext *context =
        (SequenceBenchmarkContext *)context_pointer;

    switch (context->kind) {
    case SEQUENCE_BENCHMARK_MERGE:
        return sequence_merge_sort(context->working, context->count,
                                   context->criterion,
                                   &context->statistics);
    case SEQUENCE_BENCHMARK_QUICK:
        return sequence_quick_sort(context->working, context->count,
                                   context->criterion,
                                   &context->statistics);
    case SEQUENCE_BENCHMARK_ITERATIVE_SEARCH:
        context->search_result = sequence_binary_search_iterative(
            context->source, context->count, &context->target,
            context->criterion, &context->search_index,
            &context->statistics);
        return context->search_result == SEQUENCE_SEARCH_FOUND;
    case SEQUENCE_BENCHMARK_RECURSIVE_SEARCH:
        context->search_result = sequence_binary_search_recursive(
            context->source, context->count, &context->target,
            context->criterion, &context->search_index,
            &context->statistics);
        return context->search_result == SEQUENCE_SEARCH_FOUND;
    }
    return 0;
}

static int sequence_benchmark_verify(void *context_pointer)
{
    SequenceBenchmarkContext *context =
        (SequenceBenchmarkContext *)context_pointer;
    size_t i;

    if (context->kind == SEQUENCE_BENCHMARK_ITERATIVE_SEARCH ||
        context->kind == SEQUENCE_BENCHMARK_RECURSIVE_SEARCH) {
        if (context->search_index >= context->count) {
            return 0;
        }
        if (context->criterion == SEQUENCE_SORT_BY_LENGTH) {
            return context->source[context->search_index].length ==
                   context->target.length;
        }
        return strcmp(context->source[context->search_index].identifier,
                      context->target.identifier) == 0;
    }

    for (i = 1; i < context->count; ++i) {
        if (context->criterion == SEQUENCE_SORT_BY_LENGTH) {
            if (context->working[i - 1].length >
                context->working[i].length) {
                return 0;
            }
        } else if (strcmp(context->working[i - 1].identifier,
                          context->working[i].identifier) > 0) {
            return 0;
        }
    }
    return 1;
}

static void sequence_benchmark_cleanup(void *context_pointer)
{
    SequenceBenchmarkContext *context =
        (SequenceBenchmarkContext *)context_pointer;

    free(context->working);
    context->working = NULL;
}

static int sequence_benchmark_operations(void *context_pointer,
                                         uint64_t *operations)
{
    SequenceBenchmarkContext *context =
        (SequenceBenchmarkContext *)context_pointer;

    *operations = (uint64_t)context->statistics.comparisons;
    return 1;
}

static int run_sequence_benchmarks(BenchmarkRecord *records,
                                   size_t *record_count,
                                   size_t count)
{
    DnaSequence *source = NULL;
    char (*identifiers)[SEQUENCE_BENCHMARK_IDENTIFIER_CAPACITY] = NULL;
    SequenceBenchmarkContext context;
    size_t i;
    int successful = 0;

    if (count > SIZE_MAX / sizeof(*source) ||
        count > SIZE_MAX / sizeof(*identifiers)) {
        return 0;
    }
    source = (DnaSequence *)malloc(count * sizeof(*source));
    identifiers = (char (*)[SEQUENCE_BENCHMARK_IDENTIFIER_CAPACITY])malloc(
        count * sizeof(*identifiers));
    if (source == NULL || identifiers == NULL) {
        goto cleanup;
    }
    for (i = 0; i < count; ++i) {
        (void)snprintf(identifiers[i], sizeof(identifiers[i]),
                       "seq%zu", count - i);
        source[i].identifier = identifiers[i];
        source[i].bases = NULL;
        source[i].length = (i * 37 + 11) % 211 + 1;
    }

    memset(&context, 0, sizeof(context));
    context.source = source;
    context.count = count;
    context.criterion = SEQUENCE_SORT_BY_IDENTIFIER;
    context.target = source[count / 2];
    context.kind = SEQUENCE_BENCHMARK_MERGE;
    if (!benchmark_measure("Merge Sort", count, 3,
                           sequence_benchmark_prepare,
                           sequence_benchmark_operation,
                           sequence_benchmark_verify,
                           sequence_benchmark_cleanup,
                           sequence_benchmark_operations, &context,
                           &records[*record_count])) {
        goto cleanup;
    }
    ++*record_count;

    context.kind = SEQUENCE_BENCHMARK_QUICK;
    if (!benchmark_measure("Quick Sort", count, 3,
                           sequence_benchmark_prepare,
                           sequence_benchmark_operation,
                           sequence_benchmark_verify,
                           sequence_benchmark_cleanup,
                           sequence_benchmark_operations, &context,
                           &records[*record_count])) {
        goto cleanup;
    }
    ++*record_count;

    if (!sequence_merge_sort(source, count, context.criterion, NULL)) {
        goto cleanup;
    }
    context.source = source;
    context.target = source[count / 2];
    context.kind = SEQUENCE_BENCHMARK_ITERATIVE_SEARCH;
    if (!benchmark_measure("Iterative Binary Search", count, 5,
                           sequence_benchmark_prepare,
                           sequence_benchmark_operation,
                           sequence_benchmark_verify,
                           sequence_benchmark_cleanup,
                           sequence_benchmark_operations, &context,
                           &records[*record_count])) {
        goto cleanup;
    }
    ++*record_count;

    context.kind = SEQUENCE_BENCHMARK_RECURSIVE_SEARCH;
    if (!benchmark_measure("Recursive Binary Search", count, 5,
                           sequence_benchmark_prepare,
                           sequence_benchmark_operation,
                           sequence_benchmark_verify,
                           sequence_benchmark_cleanup,
                           sequence_benchmark_operations, &context,
                           &records[*record_count])) {
        goto cleanup;
    }
    ++*record_count;
    successful = 1;

cleanup:
    free(identifiers);
    free(source);
    return successful;
}

typedef struct {
    const char *text;
    size_t text_length;
    const char *pattern;
    size_t pattern_length;
    const StringMatchResult *reference;
    StringMatchResult result;
    int algorithm;
    size_t measured_result_bytes;
    size_t result_count;
} StringBenchmarkContext;

static int string_benchmark_prepare(void *context_pointer)
{
    StringBenchmarkContext *context =
        (StringBenchmarkContext *)context_pointer;
    string_match_result_destroy(&context->result);
    return 1;
}

static int string_benchmark_operation(void *context_pointer)
{
    StringBenchmarkContext *context =
        (StringBenchmarkContext *)context_pointer;

    if (context->algorithm == 0) {
        return string_match_naive(context->text, context->text_length,
                                  context->pattern, context->pattern_length,
                                  &context->result);
    }
    if (context->algorithm == 1) {
        return string_match_rabin_karp(
            context->text, context->text_length,
            context->pattern, context->pattern_length, &context->result);
    }
    if (context->algorithm == 2) {
        return string_match_kmp(context->text, context->text_length,
                                context->pattern,
                                context->pattern_length, &context->result);
    }
    return string_match_finite_automaton(
        context->text, context->text_length, context->pattern,
        context->pattern_length, &context->result);
}

static int string_benchmark_verify(void *context_pointer)
{
    StringBenchmarkContext *context =
        (StringBenchmarkContext *)context_pointer;

    if (context->result.count != context->reference->count ||
        (context->result.count != 0 &&
         memcmp(context->result.positions, context->reference->positions,
                context->result.count *
                    sizeof(*context->result.positions)) != 0) ||
        context->result.capacity > SIZE_MAX / sizeof(size_t)) {
        return 0;
    }
    context->measured_result_bytes =
        context->result.capacity * sizeof(size_t);
    context->result_count = context->result.count;
    return 1;
}

static void string_benchmark_cleanup(void *context_pointer)
{
    StringBenchmarkContext *context =
        (StringBenchmarkContext *)context_pointer;
    string_match_result_destroy(&context->result);
}

static int run_string_benchmarks(BenchmarkRecord *records,
                                 size_t *record_count,
                                 const char *dataset,
                                 size_t text_length,
                                 size_t pattern_length,
                                 int adversarial)
{
    const char *names[] = {"Naive", "Rabin-Karp", "KMP",
                           "Finite Automaton"};
    char *text = NULL;
    char *pattern = NULL;
    StringMatchResult reference;
    StringBenchmarkContext context;
    uint32_t random_state = UINT32_C(0x71A5EED);
    size_t input_size;
    size_t i;
    int successful = 0;

    string_match_result_init(&reference);
    if (dataset == NULL || dataset[0] == '\0' || text_length == 0 ||
        pattern_length == 0 || pattern_length > text_length ||
        pattern_length > SIZE_MAX - 1 ||
        text_length > SIZE_MAX - pattern_length) {
        goto cleanup;
    }
    input_size = text_length + pattern_length;
    text = (char *)malloc(text_length);
    pattern = (char *)malloc(pattern_length + 1);
    if (text == NULL || pattern == NULL) {
        goto cleanup;
    }
    if (adversarial) {
        memset(text, 'A', text_length);
        memset(pattern, 'A', pattern_length);
        pattern[pattern_length - 1] = 'B';
    } else {
        for (i = 0; i < text_length; ++i) {
            random_state = random_state * UINT32_C(1664525) +
                           UINT32_C(1013904223);
            text[i] = "ACGT"[(random_state >> 30) & UINT32_C(3)];
        }
        memcpy(pattern, text + (text_length - pattern_length) / 2,
               pattern_length);
    }
    pattern[pattern_length] = '\0';
    if (!string_match_naive(text, text_length, pattern, pattern_length,
                            &reference)) {
        goto cleanup;
    }

    memset(&context, 0, sizeof(context));
    context.text = text;
    context.text_length = text_length;
    context.pattern = pattern;
    context.pattern_length = pattern_length;
    context.reference = &reference;
    for (i = 0; i < sizeof(names) / sizeof(names[0]); ++i) {
        context.algorithm = (int)i;
        string_match_result_init(&context.result);
        if (!benchmark_measure(names[i], input_size, 3,
                               string_benchmark_prepare,
                               string_benchmark_operation,
                               string_benchmark_verify,
                               string_benchmark_cleanup, NULL, &context,
                               &records[*record_count])) {
            string_match_result_destroy(&context.result);
            goto cleanup;
        }
        if (snprintf(records[*record_count].dataset,
                     sizeof(records[*record_count].dataset), "%s",
                     dataset) < 0 ||
            strlen(dataset) >= sizeof(records[*record_count].dataset)) {
            goto cleanup;
        }
        records[*record_count].measured_result_bytes =
            context.measured_result_bytes;
        records[*record_count].has_measured_result_bytes = 1;
        records[*record_count].result_count = context.result_count;
        records[*record_count].has_result_count = 1;
        records[*record_count].has_estimated_auxiliary_bytes = 1;
        if (context.algorithm == 2) {
            records[*record_count].estimated_auxiliary_bytes =
                pattern_length * sizeof(size_t);
        } else if (context.algorithm == 3) {
            size_t states = pattern_length + 1;
            size_t table_entries;
            size_t table_bytes;
            size_t lps_bytes;

            if (states > SIZE_MAX / 256 ||
                pattern_length > SIZE_MAX / sizeof(size_t)) {
                goto cleanup;
            }
            table_entries = states * 256;
            if (table_entries > SIZE_MAX / sizeof(size_t)) {
                goto cleanup;
            }
            table_bytes = table_entries * sizeof(size_t);
            lps_bytes = pattern_length * sizeof(size_t);
            if (lps_bytes > SIZE_MAX - table_bytes) {
                goto cleanup;
            }
            records[*record_count].estimated_auxiliary_bytes =
                table_bytes + lps_bytes;
        } else {
            records[*record_count].estimated_auxiliary_bytes = 0;
        }
        ++*record_count;
    }
    successful = 1;

cleanup:
    free(pattern);
    free(text);
    string_match_result_destroy(&reference);
    return successful;
}

typedef enum {
    KNAPSACK_BENCHMARK_DP,
    KNAPSACK_BENCHMARK_BACKTRACKING,
    KNAPSACK_BENCHMARK_FIFO,
    KNAPSACK_BENCHMARK_LIFO,
    KNAPSACK_BENCHMARK_LEAST_COST
} KnapsackBenchmarkKind;

typedef struct {
    const size_t *weights;
    const int64_t *values;
    size_t item_count;
    size_t capacity;
    int64_t expected;
    int64_t actual;
    KnapsackSearchResult result;
    KnapsackBenchmarkKind kind;
} KnapsackBenchmarkContext;

static int knapsack_benchmark_prepare(void *context_pointer)
{
    KnapsackBenchmarkContext *context =
        (KnapsackBenchmarkContext *)context_pointer;
    knapsack_search_result_destroy(&context->result);
    context->actual = 0;
    return 1;
}

static int knapsack_benchmark_operation(void *context_pointer)
{
    KnapsackBenchmarkContext *context =
        (KnapsackBenchmarkContext *)context_pointer;

    switch (context->kind) {
    case KNAPSACK_BENCHMARK_DP:
        return dp_knapsack_01(context->weights, context->values,
                              context->item_count, context->capacity,
                              &context->actual);
    case KNAPSACK_BENCHMARK_BACKTRACKING:
        return backtracking_knapsack_01(
            context->weights, context->values, context->item_count,
            context->capacity, &context->result);
    case KNAPSACK_BENCHMARK_FIFO:
        return branch_bound_fifo_knapsack_01(
            context->weights, context->values, context->item_count,
            context->capacity, &context->result);
    case KNAPSACK_BENCHMARK_LIFO:
        return branch_bound_lifo_knapsack_01(
            context->weights, context->values, context->item_count,
            context->capacity, &context->result);
    case KNAPSACK_BENCHMARK_LEAST_COST:
        return branch_bound_least_cost_knapsack_01(
            context->weights, context->values, context->item_count,
            context->capacity, &context->result);
    }
    return 0;
}

static int knapsack_benchmark_verify(void *context_pointer)
{
    KnapsackBenchmarkContext *context =
        (KnapsackBenchmarkContext *)context_pointer;
    int64_t value = context->kind == KNAPSACK_BENCHMARK_DP
                        ? context->actual
                        : context->result.maximum_value;
    return value == context->expected;
}

static void knapsack_benchmark_cleanup(void *context_pointer)
{
    KnapsackBenchmarkContext *context =
        (KnapsackBenchmarkContext *)context_pointer;
    knapsack_search_result_destroy(&context->result);
}

static int run_knapsack_benchmarks(BenchmarkRecord *records,
                                   size_t *record_count,
                                   size_t item_count)
{
    const char *names[] = {"DP Knapsack", "Backtracking Knapsack",
                           "FIFO Branch and Bound",
                           "LIFO Branch and Bound",
                           "Least-Cost Branch and Bound"};
    size_t *weights = NULL;
    int64_t *values = NULL;
    size_t total_weight = 0;
    size_t capacity;
    size_t i;
    int64_t expected;
    int successful = 0;
    KnapsackBenchmarkContext context = {0};

    if (item_count == 0 ||
        item_count > SIZE_MAX / sizeof(*weights) ||
        item_count > SIZE_MAX / sizeof(*values)) {
        return 0;
    }
    weights = (size_t *)malloc(item_count * sizeof(*weights));
    values = (int64_t *)malloc(item_count * sizeof(*values));
    if (weights == NULL || values == NULL) {
        goto cleanup;
    }
    for (i = 0; i < item_count; ++i) {
        weights[i] = 2 + (i * 13 + 7) % 19;
        values[i] = (int64_t)(5 + (i * 17 + 3) % 41);
        if (weights[i] > SIZE_MAX - total_weight) {
            goto cleanup;
        }
        total_weight += weights[i];
    }
    capacity = total_weight / 3;
    if (!dp_knapsack_01(weights, values, item_count, capacity, &expected)) {
        goto cleanup;
    }

    memset(&context, 0, sizeof(context));
    context.weights = weights;
    context.values = values;
    context.item_count = item_count;
    context.capacity = capacity;
    context.expected = expected;
    for (i = 0; i < 5; ++i) {
        context.kind = (KnapsackBenchmarkKind)i;
        if (!benchmark_measure(names[i], item_count, 1,
                               knapsack_benchmark_prepare,
                               knapsack_benchmark_operation,
                               knapsack_benchmark_verify,
                               knapsack_benchmark_cleanup, NULL, &context,
                               &records[*record_count])) {
            goto cleanup;
        }
        records[*record_count].solution_value = expected;
        records[*record_count].has_solution_value = 1;
        ++*record_count;
    }
    successful = 1;

cleanup:
    knapsack_search_result_destroy(&context.result);
    free(values);
    free(weights);
    return successful;
}

typedef enum {
    GRAPH_BENCHMARK_DIJKSTRA,
    GRAPH_BENCHMARK_PRIM,
    GRAPH_BENCHMARK_KRUSKAL
} GraphBenchmarkKind;

typedef struct {
    WeightedGraph graph;
    WeightedEdge *edges;
    size_t edge_count;
    GraphBenchmarkKind kind;
    ShortestPathResult shortest;
    PrimResult prim;
    KruskalResult kruskal;
    double expected_mst;
} GraphBenchmarkContext;

static int graph_benchmark_prepare(void *context_pointer)
{
    GraphBenchmarkContext *context =
        (GraphBenchmarkContext *)context_pointer;
    greedy_shortest_path_destroy(&context->shortest);
    greedy_prim_destroy(&context->prim);
    greedy_kruskal_destroy(&context->kruskal);
    return 1;
}

static int graph_benchmark_operation(void *context_pointer)
{
    GraphBenchmarkContext *context =
        (GraphBenchmarkContext *)context_pointer;

    if (context->kind == GRAPH_BENCHMARK_DIJKSTRA) {
        return greedy_dijkstra_matrix(&context->graph, 0,
                                      &context->shortest);
    }
    if (context->kind == GRAPH_BENCHMARK_PRIM) {
        return greedy_prim_matrix(&context->graph, 0, &context->prim);
    }
    return greedy_kruskal(context->graph.vertex_count,
                          context->edges, context->edge_count,
                          &context->kruskal);
}

static int graph_benchmark_verify(void *context_pointer)
{
    GraphBenchmarkContext *context =
        (GraphBenchmarkContext *)context_pointer;
    size_t i;

    if (context->kind == GRAPH_BENCHMARK_DIJKSTRA) {
        if (context->shortest.vertex_count != context->graph.vertex_count ||
            context->shortest.distances[0] != 0.0) {
            return 0;
        }
        for (i = 0; i < context->graph.vertex_count; ++i) {
            if (!isfinite(context->shortest.distances[i])) {
                return 0;
            }
        }
        return 1;
    }
    if (context->kind == GRAPH_BENCHMARK_PRIM) {
        return context->prim.edge_count + 1 ==
                   context->graph.vertex_count &&
               context->prim.total_weight == context->expected_mst;
    }
    return context->kruskal.connected &&
           context->kruskal.total_weight == context->expected_mst;
}

static void graph_benchmark_cleanup(void *context_pointer)
{
    GraphBenchmarkContext *context =
        (GraphBenchmarkContext *)context_pointer;
    greedy_shortest_path_destroy(&context->shortest);
    greedy_prim_destroy(&context->prim);
    greedy_kruskal_destroy(&context->kruskal);
}

static int run_graph_benchmarks(BenchmarkRecord *records,
                                size_t *record_count,
                                size_t vertex_count)
{
    const char *names[] = {"Dijkstra", "Prim", "Kruskal"};
    GraphBenchmarkContext context;
    size_t matrix_count;
    size_t edge_count;
    size_t i;
    size_t j;
    size_t edge_index = 0;
    int successful = 0;

    if (vertex_count < 2 ||
        vertex_count > SIZE_MAX / vertex_count) {
        return 0;
    }
    matrix_count = vertex_count * vertex_count;
    if (matrix_count > SIZE_MAX / sizeof(double) ||
        vertex_count > SIZE_MAX / (vertex_count - 1)) {
        return 0;
    }
    edge_count = vertex_count * (vertex_count - 1) / 2;
    if (edge_count > SIZE_MAX / sizeof(*context.edges)) {
        return 0;
    }
    memset(&context, 0, sizeof(context));
    context.graph.vertex_count = vertex_count;
    context.graph.adjacency = (double *)malloc(matrix_count * sizeof(double));
    context.edges = (WeightedEdge *)malloc(edge_count * sizeof(*context.edges));
    context.edge_count = edge_count;
    if (context.graph.adjacency == NULL || context.edges == NULL) {
        goto cleanup;
    }

    for (i = 0; i < vertex_count; ++i) {
        for (j = 0; j < vertex_count; ++j) {
            if (i == j) {
                context.graph.adjacency[i * vertex_count + j] = 0.0;
            } else {
                size_t low = i < j ? i : j;
                size_t high = i < j ? j : i;
                double weight = (double)(1 + (low * 31 + high * 17) % 97);

                context.graph.adjacency[i * vertex_count + j] = weight;
                if (i < j) {
                    context.edges[edge_index].from = i;
                    context.edges[edge_index].to = j;
                    context.edges[edge_index].weight = weight;
                    ++edge_index;
                }
            }
        }
    }
    if (!greedy_prim_matrix(&context.graph, 0, &context.prim)) {
        goto cleanup;
    }
    context.expected_mst = context.prim.total_weight;
    greedy_prim_destroy(&context.prim);

    for (i = 0; i < 3; ++i) {
        context.kind = (GraphBenchmarkKind)i;
        if (!benchmark_measure(names[i], vertex_count, 3,
                               graph_benchmark_prepare,
                               graph_benchmark_operation,
                               graph_benchmark_verify,
                               graph_benchmark_cleanup, NULL, &context,
                               &records[*record_count])) {
            goto cleanup;
        }
        ++*record_count;
    }
    successful = 1;

cleanup:
    graph_benchmark_cleanup(&context);
    free(context.edges);
    free(context.graph.adjacency);
    return successful;
}

typedef struct {
    FordFulkersonNetwork network;
    FordFulkersonEdge *edges;
    FordFulkersonResult result;
} FlowBenchmarkContext;

static int flow_benchmark_prepare(void *context_pointer)
{
    FlowBenchmarkContext *context =
        (FlowBenchmarkContext *)context_pointer;
    ford_fulkerson_result_destroy(&context->result);
    return 1;
}

static int flow_benchmark_operation(void *context_pointer)
{
    FlowBenchmarkContext *context =
        (FlowBenchmarkContext *)context_pointer;
    return ford_fulkerson_max_flow(&context->network, 0,
                                   context->network.vertex_count - 1,
                                   &context->result);
}

static int flow_benchmark_verify(void *context_pointer)
{
    FlowBenchmarkContext *context =
        (FlowBenchmarkContext *)context_pointer;
    return context->result.max_flow == context->result.cut_capacity;
}

static void flow_benchmark_cleanup(void *context_pointer)
{
    FlowBenchmarkContext *context =
        (FlowBenchmarkContext *)context_pointer;
    ford_fulkerson_result_destroy(&context->result);
}

static int run_flow_benchmark(BenchmarkRecord *records,
                              size_t *record_count,
                              size_t vertex_count)
{
    FlowBenchmarkContext context;
    size_t edge_count;
    size_t i;
    size_t edge_index = 0;
    int successful = 0;

    if (vertex_count < 2 || vertex_count > SIZE_MAX / 2) {
        return 0;
    }
    edge_count = 2 * vertex_count - 3;
    if (edge_count > SIZE_MAX / sizeof(*context.edges)) {
        return 0;
    }
    memset(&context, 0, sizeof(context));
    context.edges = (FordFulkersonEdge *)malloc(
        edge_count * sizeof(*context.edges));
    if (context.edges == NULL) {
        return 0;
    }
    for (i = 0; i + 1 < vertex_count; ++i) {
        context.edges[edge_index].from = i;
        context.edges[edge_index].to = i + 1;
        context.edges[edge_index].capacity = (int64_t)(10 + (i % 7));
        ++edge_index;
    }
    for (i = 0; i + 2 < vertex_count; ++i) {
        context.edges[edge_index].from = i;
        context.edges[edge_index].to = i + 2;
        context.edges[edge_index].capacity = (int64_t)(5 + (i % 5));
        ++edge_index;
    }
    context.network.vertex_count = vertex_count;
    context.network.edge_count = edge_count;
    context.network.edges = context.edges;

    successful = benchmark_measure("Ford-Fulkerson", vertex_count, 3,
                                   flow_benchmark_prepare,
                                   flow_benchmark_operation,
                                   flow_benchmark_verify,
                                   flow_benchmark_cleanup, NULL, &context,
                                   &records[*record_count]);
    if (successful) {
        ++*record_count;
    }
    flow_benchmark_cleanup(&context);
    free(context.edges);
    return successful;
}

int benchmark_run_suite(BenchmarkRecord *records,
                        size_t record_capacity,
                        size_t *record_count)
{
    static const size_t sequence_sizes[] = {64, 256};
    static const struct {
        const char *name;
        size_t text_length;
        size_t pattern_length;
        int adversarial;
    } string_datasets[] = {
        {"small", 256, 8, 0},
        {"medium", 8192, 24, 0},
        {"large", 65536, 32, 0},
        {"adversarial", 32768, 64, 1}
    };
    static const size_t knapsack_sizes[] = {10, 14};
    static const size_t graph_sizes[] = {16, 32};
    size_t local_count = 0;
    size_t i;

    if (records == NULL || record_count == NULL ||
        record_capacity < BENCHMARK_SUITE_MAX_RECORDS) {
        return 0;
    }
    *record_count = 0;
    for (i = 0; i < sizeof(sequence_sizes) / sizeof(sequence_sizes[0]); ++i) {
        if (!run_sequence_benchmarks(records, &local_count,
                                     sequence_sizes[i])) {
            goto failure;
        }
    }
    for (i = 0; i < sizeof(string_datasets) / sizeof(string_datasets[0]); ++i) {
        if (!run_string_benchmarks(records, &local_count,
                                  string_datasets[i].name,
                                  string_datasets[i].text_length,
                                  string_datasets[i].pattern_length,
                                  string_datasets[i].adversarial)) {
            goto failure;
        }
    }
    for (i = 0; i < sizeof(knapsack_sizes) / sizeof(knapsack_sizes[0]); ++i) {
        if (!run_knapsack_benchmarks(records, &local_count,
                                     knapsack_sizes[i])) {
            goto failure;
        }
    }
    for (i = 0; i < sizeof(graph_sizes) / sizeof(graph_sizes[0]); ++i) {
        if (!run_graph_benchmarks(records, &local_count, graph_sizes[i]) ||
            !run_flow_benchmark(records, &local_count, graph_sizes[i])) {
            goto failure;
        }
    }
    *record_count = local_count;
    return 1;

failure:
    *record_count = local_count;
    return 0;
}
