#ifndef GENOME_ALGORITHMS_BACKTRACKING_ALGORITHMS_H
#define GENOME_ALGORITHMS_BACKTRACKING_ALGORITHMS_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t *indices;
    size_t index_count;
} SumOfSubsetsSubset;

typedef struct {
    SumOfSubsetsSubset *subsets;
    size_t subset_count;
} SumOfSubsetsResult;

int backtracking_sum_of_subsets(const int64_t *values,
                                size_t item_count,
                                int64_t target,
                                SumOfSubsetsResult *result);
void backtracking_sum_of_subsets_result_destroy(SumOfSubsetsResult *result);

typedef struct {
    int64_t maximum_value;
    size_t *selected_indices;
    size_t selected_count;
} KnapsackSearchResult;

int backtracking_knapsack_01(const size_t *weights,
                             const int64_t *values,
                             size_t item_count,
                             size_t capacity,
                             KnapsackSearchResult *result);

int branch_bound_fifo_knapsack_01(const size_t *weights,
                                  const int64_t *values,
                                  size_t item_count,
                                  size_t capacity,
                                  KnapsackSearchResult *result);
int branch_bound_lifo_knapsack_01(const size_t *weights,
                                  const int64_t *values,
                                  size_t item_count,
                                  size_t capacity,
                                  KnapsackSearchResult *result);
int branch_bound_least_cost_knapsack_01(const size_t *weights,
                                        const int64_t *values,
                                        size_t item_count,
                                        size_t capacity,
                                        KnapsackSearchResult *result);
void knapsack_search_result_destroy(KnapsackSearchResult *result);

#endif
