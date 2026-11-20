#include "genome/algorithms/backtracking_algorithms.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int64_t value;
    size_t original_index;
} SubsetItem;

typedef struct {
    SubsetItem *items;
    size_t item_count;
    size_t *current_indices;
    size_t current_count;
    SumOfSubsetsResult *result;
    size_t result_capacity;
    int failed;
} SumOfSubsetsSearch;

static int compare_subset_items(const void *left_value,
                                const void *right_value)
{
    const SubsetItem *left = (const SubsetItem *)left_value;
    const SubsetItem *right = (const SubsetItem *)right_value;

    if (left->value < right->value) {
        return -1;
    }
    if (left->value > right->value) {
        return 1;
    }
    if (left->original_index < right->original_index) {
        return -1;
    }
    return left->original_index > right->original_index;
}

static int compare_size_values(const void *left_value,
                               const void *right_value)
{
    size_t left = *(const size_t *)left_value;
    size_t right = *(const size_t *)right_value;

    return left < right ? -1 : left > right;
}

static int sum_of_subsets_reserve(SumOfSubsetsSearch *search)
{
    size_t new_capacity;
    SumOfSubsetsSubset *new_subsets;

    if (search->result->subset_count == SIZE_MAX) {
        return 0;
    }
    if (search->result->subset_count < search->result_capacity) {
        return 1;
    }
    if (search->result_capacity == 0) {
        new_capacity = 4;
    } else {
        if (search->result_capacity > SIZE_MAX / 2) {
            return 0;
        }
        new_capacity = search->result_capacity * 2;
    }
    if (new_capacity > SIZE_MAX / sizeof(*new_subsets)) {
        return 0;
    }
    new_subsets = (SumOfSubsetsSubset *)realloc(
        search->result->subsets, new_capacity * sizeof(*new_subsets));
    if (new_subsets == NULL) {
        return 0;
    }
    search->result->subsets = new_subsets;
    search->result_capacity = new_capacity;
    return 1;
}

static int sum_of_subsets_record(SumOfSubsetsSearch *search)
{
    SumOfSubsetsSubset subset = {NULL, search->current_count};

    if (!sum_of_subsets_reserve(search)) {
        return 0;
    }
    if (subset.index_count != 0) {
        if (subset.index_count > SIZE_MAX / sizeof(*subset.indices)) {
            return 0;
        }
        subset.indices = (size_t *)malloc(
            subset.index_count * sizeof(*subset.indices));
        if (subset.indices == NULL) {
            return 0;
        }
        memcpy(subset.indices, search->current_indices,
               subset.index_count * sizeof(*subset.indices));
        qsort(subset.indices, subset.index_count,
              sizeof(*subset.indices), compare_size_values);
    }
    search->result->subsets[search->result->subset_count++] = subset;
    return 1;
}

static void sum_of_subsets_search(SumOfSubsetsSearch *search,
                                  size_t start,
                                  int64_t remaining)
{
    size_t i;

    if (search->failed) {
        return;
    }
    if (remaining == 0) {
        if (!sum_of_subsets_record(search)) {
            search->failed = 1;
        }
        return;
    }

    for (i = start; i < search->item_count; ++i) {
        int64_t value = search->items[i].value;

        if (value > remaining) {
            break;
        }
        search->current_indices[search->current_count++] =
            search->items[i].original_index;
        sum_of_subsets_search(search, i + 1, remaining - value);
        --search->current_count;
        if (search->failed) {
            return;
        }
    }
}

void backtracking_sum_of_subsets_result_destroy(SumOfSubsetsResult *result)
{
    size_t i;

    if (result == NULL) {
        return;
    }
    for (i = 0; i < result->subset_count; ++i) {
        free(result->subsets[i].indices);
    }
    free(result->subsets);
    result->subsets = NULL;
    result->subset_count = 0;
}

int backtracking_sum_of_subsets(const int64_t *values,
                                size_t item_count,
                                int64_t target,
                                SumOfSubsetsResult *result)
{
    SumOfSubsetsSearch search;
    SumOfSubsetsResult completed = {NULL, 0};
    size_t i;

    if (result == NULL || target < 0 ||
        (item_count != 0 && values == NULL) ||
        item_count > SIZE_MAX / sizeof(*search.items) ||
        item_count > SIZE_MAX / sizeof(*search.current_indices)) {
        return 0;
    }

    memset(&search, 0, sizeof(search));
    search.item_count = item_count;
    search.result = &completed;
    if (item_count != 0) {
        search.items =
            (SubsetItem *)malloc(item_count * sizeof(*search.items));
        search.current_indices = (size_t *)malloc(
            item_count * sizeof(*search.current_indices));
        if (search.items == NULL || search.current_indices == NULL) {
            search.failed = 1;
        } else {
            for (i = 0; i < item_count; ++i) {
                if (values[i] <= 0) {
                    search.failed = 1;
                    break;
                }
                search.items[i].value = values[i];
                search.items[i].original_index = i;
            }
            if (!search.failed) {
                qsort(search.items, item_count, sizeof(*search.items),
                      compare_subset_items);
            }
        }
    }

    if (!search.failed) {
        if (target == 0) {
            if (!sum_of_subsets_record(&search)) {
                search.failed = 1;
            }
        } else {
            sum_of_subsets_search(&search, 0, target);
        }
    }
    free(search.current_indices);
    free(search.items);
    if (search.failed) {
        backtracking_sum_of_subsets_result_destroy(&completed);
        return 0;
    }

    backtracking_sum_of_subsets_result_destroy(result);
    *result = completed;
    return 1;
}

typedef struct {
    const size_t *weights;
    const int64_t *values;
    size_t item_count;
    size_t capacity;
    size_t *current_indices;
    size_t *best_indices;
    size_t current_count;
    size_t best_count;
    int64_t best_value;
    int failed;
} BacktrackingKnapsackSearch;

static void backtracking_knapsack_visit(BacktrackingKnapsackSearch *search,
                                        size_t item_index,
                                        size_t current_weight,
                                        int64_t current_value)
{
    int64_t item_value;

    if (search->failed) {
        return;
    }
    if (item_index == search->item_count) {
        if (current_value > search->best_value) {
            search->best_value = current_value;
            search->best_count = search->current_count;
            if (search->best_count != 0) {
                memcpy(search->best_indices, search->current_indices,
                       search->best_count * sizeof(*search->best_indices));
            }
        }
        return;
    }

    backtracking_knapsack_visit(search, item_index + 1,
                                current_weight, current_value);
    if (search->failed ||
        search->weights[item_index] >
            search->capacity - current_weight) {
        return;
    }

    item_value = search->values[item_index];
    if (item_value > INT64_MAX - current_value) {
        search->failed = 1;
        return;
    }
    search->current_indices[search->current_count++] = item_index;
    backtracking_knapsack_visit(search, item_index + 1,
                                current_weight + search->weights[item_index],
                                current_value + item_value);
    --search->current_count;
}

void knapsack_search_result_destroy(KnapsackSearchResult *result)
{
    if (result == NULL) {
        return;
    }
    free(result->selected_indices);
    result->maximum_value = 0;
    result->selected_indices = NULL;
    result->selected_count = 0;
}

static int knapsack_inputs_are_valid(const size_t *weights,
                                    const int64_t *values,
                                    size_t item_count)
{
    size_t i;

    if (item_count > SIZE_MAX / sizeof(*weights) ||
        item_count > SIZE_MAX / sizeof(*values) ||
        (item_count != 0 && (weights == NULL || values == NULL))) {
        return 0;
    }
    for (i = 0; i < item_count; ++i) {
        if (values[i] < 0) {
            return 0;
        }
    }
    return 1;
}

static int knapsack_result_replace(KnapsackSearchResult *result,
                                   int64_t maximum_value,
                                   const size_t *indices,
                                   size_t index_count)
{
    KnapsackSearchResult completed = {maximum_value, NULL, index_count};

    if (index_count != 0) {
        if (index_count > SIZE_MAX / sizeof(*completed.selected_indices)) {
            return 0;
        }
        completed.selected_indices = (size_t *)malloc(
            index_count * sizeof(*completed.selected_indices));
        if (completed.selected_indices == NULL) {
            return 0;
        }
        memcpy(completed.selected_indices, indices,
               index_count * sizeof(*completed.selected_indices));
        qsort(completed.selected_indices, index_count,
              sizeof(*completed.selected_indices), compare_size_values);
    }
    knapsack_search_result_destroy(result);
    *result = completed;
    return 1;
}

int backtracking_knapsack_01(const size_t *weights,
                             const int64_t *values,
                             size_t item_count,
                             size_t capacity,
                             KnapsackSearchResult *result)
{
    BacktrackingKnapsackSearch search;
    int success;

    if (result == NULL ||
        item_count > SIZE_MAX / sizeof(*search.current_indices) ||
        !knapsack_inputs_are_valid(weights, values, item_count)) {
        return 0;
    }
    if (item_count == 0) {
        return knapsack_result_replace(result, 0, NULL, 0);
    }
    memset(&search, 0, sizeof(search));
    search.weights = weights;
    search.values = values;
    search.item_count = item_count;
    search.capacity = capacity;
    search.current_indices =
        (size_t *)malloc(item_count * sizeof(*search.current_indices));
    search.best_indices =
        (size_t *)malloc(item_count * sizeof(*search.best_indices));
    if (search.current_indices == NULL || search.best_indices == NULL) {
        free(search.best_indices);
        free(search.current_indices);
        return 0;
    }

    backtracking_knapsack_visit(&search, 0, 0, 0);
    success = !search.failed &&
              knapsack_result_replace(result, search.best_value,
                                      search.best_indices,
                                      search.best_count);
    free(search.best_indices);
    free(search.current_indices);
    return success;
}

typedef struct {
    size_t original_index;
    size_t weight;
    int64_t value;
} BoundItem;

typedef struct {
    size_t level;
    size_t weight;
    int64_t value;
    int64_t bound;
    size_t parent;
    size_t original_index;
    int included;
} BoundNode;

typedef enum {
    OPEN_FIFO,
    OPEN_LIFO,
    OPEN_BEST_BOUND
} OpenPolicy;

typedef struct {
    size_t *indices;
    size_t count;
    size_t capacity;
    size_t head;
    OpenPolicy policy;
} OpenSet;

static int compare_positive_fractions(uint64_t left_numerator,
                                      uint64_t left_denominator,
                                      uint64_t right_numerator,
                                      uint64_t right_denominator)
{
    int direction = 1;

    for (;;) {
        uint64_t left_whole =
            left_numerator / left_denominator;
        uint64_t right_whole =
            right_numerator / right_denominator;
        uint64_t left_remainder =
            left_numerator % left_denominator;
        uint64_t right_remainder =
            right_numerator % right_denominator;

        if (left_whole != right_whole) {
            return direction *
                   (left_whole < right_whole ? -1 : 1);
        }
        if (left_remainder == 0 || right_remainder == 0) {
            if (left_remainder == 0 && right_remainder == 0) {
                return 0;
            }
            return direction *
                   (left_remainder == 0 ? -1 : 1);
        }

        left_numerator = left_denominator;
        left_denominator = left_remainder;
        right_numerator = right_denominator;
        right_denominator = right_remainder;
        direction = -direction;
    }
}

static int compare_bound_items(const void *left_value,
                               const void *right_value)
{
    const BoundItem *left = (const BoundItem *)left_value;
    const BoundItem *right = (const BoundItem *)right_value;
    int ratio_comparison;

    if (left->weight == 0 || right->weight == 0) {
        int left_infinite = left->weight == 0 && left->value > 0;
        int right_infinite = right->weight == 0 && right->value > 0;

        if (left_infinite != right_infinite) {
            return left_infinite ? -1 : 1;
        }
        if (left->weight == 0 && right->weight != 0 &&
            left->value == 0) {
            ratio_comparison =
                right->value == 0 ? 0 : -1;
        } else if (right->weight == 0 && left->weight != 0 &&
                   right->value == 0) {
            ratio_comparison =
                left->value == 0 ? 0 : 1;
        } else {
            ratio_comparison = 0;
        }
    } else {
        ratio_comparison = compare_positive_fractions(
            (uint64_t)left->value, (uint64_t)left->weight,
            (uint64_t)right->value, (uint64_t)right->weight);
    }
    if (ratio_comparison != 0) {
        return ratio_comparison > 0 ? -1 : 1;
    }
    if (left->original_index < right->original_index) {
        return -1;
    }
    return left->original_index > right->original_index;
}

static int prepare_bound_items(const size_t *weights,
                               const int64_t *values,
                               size_t item_count,
                               BoundItem **items)
{
    BoundItem *ordered;
    size_t i;

    *items = NULL;
    if (item_count > SIZE_MAX / sizeof(*ordered) ||
        !knapsack_inputs_are_valid(weights, values, item_count)) {
        return 0;
    }
    if (item_count == 0) {
        return 1;
    }

    ordered = (BoundItem *)malloc(item_count * sizeof(*ordered));
    if (ordered == NULL) {
        return 0;
    }
    for (i = 0; i < item_count; ++i) {
        ordered[i].original_index = i;
        ordered[i].weight = weights[i];
        ordered[i].value = values[i];
    }
    qsort(ordered, item_count, sizeof(*ordered), compare_bound_items);
    *items = ordered;
    return 1;
}

static uint64_t multiply_divide_ceil(uint64_t value,
                                    uint64_t multiplier,
                                    uint64_t divisor)
{
    uint64_t quotient = 0;
    uint64_t remainder = 0;
    int bit;

    for (bit = 63; bit >= 0; --bit) {
        int carry;

        quotient *= 2;
        if (remainder >= divisor - remainder) {
            remainder -= divisor - remainder;
            carry = 1;
        } else {
            remainder += remainder;
            carry = 0;
        }
        quotient += (uint64_t)carry;

        if (((value >> bit) & UINT64_C(1)) != 0) {
            if (remainder >= divisor - multiplier) {
                remainder -= divisor - multiplier;
                ++quotient;
            } else {
                remainder += multiplier;
            }
        }
    }
    if (remainder != 0) {
        ++quotient;
    }
    return quotient;
}

static int64_t add_bound_saturated(int64_t current, uint64_t addition)
{
    uint64_t limit = (uint64_t)INT64_MAX;
    uint64_t base = (uint64_t)current;

    if (addition >= limit - base) {
        return INT64_MAX;
    }
    return (int64_t)(base + addition);
}

static int64_t density_upper_bound(const BoundItem *items,
                                   const int64_t *suffix_values,
                                   const int64_t *suffix_zero_weight_values,
                                   const size_t *maximum_density_item,
                                   size_t level,
                                   size_t current_weight,
                                   int64_t current_value,
                                   size_t capacity)
{
    uint64_t remaining = (uint64_t)(capacity - current_weight);
    uint64_t density_bound =
        (uint64_t)suffix_zero_weight_values[level];
    int64_t maximum_additional;
    size_t density_item = maximum_density_item[level];

    if (density_item != SIZE_MAX && remaining != 0) {
        uint64_t value = (uint64_t)items[density_item].value;
        uint64_t weight = (uint64_t)items[density_item].weight;
        uint64_t whole = remaining / weight;
        uint64_t remainder = remaining % weight;
        uint64_t limit = (uint64_t)INT64_MAX;
        uint64_t ratio_bound;

        if (whole != 0 && value > limit / whole) {
            ratio_bound = limit;
        } else {
            ratio_bound = value * whole;
            if (remainder != 0) {
                uint64_t fractional = multiply_divide_ceil(
                    value, remainder, weight);

                if (fractional >= limit - ratio_bound) {
                    ratio_bound = limit;
                } else {
                    ratio_bound += fractional;
                }
            }
        }
        if (ratio_bound >= limit - density_bound) {
            density_bound = limit;
        } else {
            density_bound += ratio_bound;
        }
    }

    maximum_additional =
        suffix_values[level] < (int64_t)density_bound
            ? suffix_values[level]
            : (int64_t)density_bound;
    return add_bound_saturated(current_value,
                               (uint64_t)maximum_additional);
}

static int node_has_higher_priority(const BoundNode *nodes,
                                   size_t left_index,
                                   size_t right_index)
{
    const BoundNode *left = &nodes[left_index];
    const BoundNode *right = &nodes[right_index];

    if (left->bound != right->bound) {
        return left->bound > right->bound;
    }
    if (left->value != right->value) {
        return left->value > right->value;
    }
    return left_index < right_index;
}

static int bound_cannot_improve(int64_t bound, int64_t best_value)
{
    return bound < best_value ||
           (bound == best_value && best_value < INT64_MAX);
}

static int open_set_reserve(OpenSet *open_set)
{
    size_t new_capacity;
    size_t *new_indices;

    if (open_set->count < open_set->capacity) {
        return 1;
    }
    if (open_set->capacity == 0) {
        new_capacity = 16;
    } else {
        if (open_set->capacity > SIZE_MAX / 2) {
            return 0;
        }
        new_capacity = open_set->capacity * 2;
    }
    if (new_capacity > SIZE_MAX / sizeof(*new_indices)) {
        return 0;
    }
    new_indices = (size_t *)realloc(
        open_set->indices, new_capacity * sizeof(*new_indices));
    if (new_indices == NULL) {
        return 0;
    }
    open_set->indices = new_indices;
    open_set->capacity = new_capacity;
    return 1;
}

static int open_set_push(OpenSet *open_set,
                         const BoundNode *nodes,
                         size_t node_index)
{
    size_t position;

    if (!open_set_reserve(open_set)) {
        return 0;
    }
    position = open_set->count++;
    open_set->indices[position] = node_index;
    if (open_set->policy == OPEN_BEST_BOUND) {
        while (position != 0) {
            size_t parent = (position - 1) / 2;
            size_t temporary;

            if (!node_has_higher_priority(
                    nodes, open_set->indices[position],
                    open_set->indices[parent])) {
                break;
            }
            temporary = open_set->indices[position];
            open_set->indices[position] = open_set->indices[parent];
            open_set->indices[parent] = temporary;
            position = parent;
        }
    }
    return 1;
}

static int open_set_is_empty(const OpenSet *open_set)
{
    if (open_set->policy == OPEN_FIFO) {
        return open_set->head == open_set->count;
    }
    return open_set->count == 0;
}

static size_t open_set_pop(OpenSet *open_set, const BoundNode *nodes)
{
    size_t answer;

    if (open_set->policy == OPEN_FIFO) {
        return open_set->indices[open_set->head++];
    }
    if (open_set->policy == OPEN_LIFO) {
        return open_set->indices[--open_set->count];
    }

    answer = open_set->indices[0];
    --open_set->count;
    if (open_set->count != 0) {
        size_t position = 0;
        size_t last = open_set->indices[open_set->count];

        for (;;) {
            size_t left = position * 2 + 1;
            size_t right = left + 1;
            size_t best = position;

            if (left < open_set->count &&
                node_has_higher_priority(
                    nodes, open_set->indices[left],
                    best == position ? last : open_set->indices[best])) {
                best = left;
            }
            if (right < open_set->count &&
                node_has_higher_priority(
                    nodes, open_set->indices[right],
                    best == position ? last : open_set->indices[best])) {
                best = right;
            }
            if (best == position) {
                break;
            }
            open_set->indices[position] = open_set->indices[best];
            position = best;
        }
        open_set->indices[position] = last;
    }
    return answer;
}

static int bound_nodes_reserve(BoundNode **nodes,
                               size_t *capacity,
                               size_t required)
{
    size_t new_capacity;
    BoundNode *new_nodes;

    if (required <= *capacity) {
        return 1;
    }
    new_capacity = *capacity == 0 ? 32 : *capacity;
    while (new_capacity < required) {
        if (new_capacity > SIZE_MAX / 2) {
            new_capacity = required;
            break;
        }
        new_capacity *= 2;
    }
    if (new_capacity > SIZE_MAX / sizeof(*new_nodes)) {
        return 0;
    }
    new_nodes = (BoundNode *)realloc(
        *nodes, new_capacity * sizeof(*new_nodes));
    if (new_nodes == NULL) {
        return 0;
    }
    *nodes = new_nodes;
    *capacity = new_capacity;
    return 1;
}

static int bound_node_append(BoundNode **nodes,
                             size_t *node_count,
                             size_t *node_capacity,
                             const BoundNode *node,
                             size_t *node_index)
{
    if (*node_count == SIZE_MAX ||
        !bound_nodes_reserve(nodes, node_capacity, *node_count + 1)) {
        return 0;
    }
    *node_index = *node_count;
    (*nodes)[(*node_count)++] = *node;
    return 1;
}

static int branch_bound_add_child(const BoundItem *items,
                                  const int64_t *suffix_values,
                                  const int64_t *suffix_zero_weight_values,
                                  const size_t *maximum_density_item,
                                  size_t capacity,
                                  const BoundNode *parent,
                                  size_t parent_index,
                                  size_t item_index,
                                  int include,
                                  int64_t best_value,
                                  BoundNode **nodes,
                                  size_t *node_count,
                                  size_t *node_capacity,
                                  OpenSet *open_set)
{
    BoundNode child;
    size_t child_index;

    child.level = parent->level + 1;
    child.weight = parent->weight;
    child.value = parent->value;
    child.parent = parent_index;
    child.original_index = items[item_index].original_index;
    child.included = include;

    if (include) {
        if (items[item_index].weight > capacity - parent->weight) {
            return 1;
        }
        if (items[item_index].value > INT64_MAX - parent->value) {
            return 0;
        }
        child.weight += items[item_index].weight;
        child.value += items[item_index].value;
    }
    child.bound = density_upper_bound(
        items, suffix_values, suffix_zero_weight_values,
        maximum_density_item, child.level, child.weight, child.value,
        capacity);
    if (bound_cannot_improve(child.bound, best_value)) {
        return 1;
    }
    if (!bound_node_append(nodes, node_count, node_capacity, &child,
                           &child_index) ||
        !open_set_push(open_set, *nodes, child_index)) {
        return 0;
    }
    return 1;
}

static int branch_bound_knapsack(const size_t *weights,
                                 const int64_t *values,
                                 size_t item_count,
                                 size_t capacity,
                                 OpenPolicy policy,
                                 KnapsackSearchResult *result)
{
    BoundItem *items = NULL;
    int64_t *suffix_values = NULL;
    int64_t *suffix_zero_weight_values = NULL;
    size_t *maximum_density_item = NULL;
    BoundNode *nodes = NULL;
    size_t node_count = 0;
    size_t node_capacity = 0;
    size_t best_node = 0;
    size_t *selected_indices = NULL;
    size_t selected_count = 0;
    OpenSet open_set = {NULL, 0, 0, 0, policy};
    BoundNode root;
    int64_t best_value = 0;
    size_t i;
    int success = 0;

    if (result == NULL ||
        !prepare_bound_items(weights, values, item_count, &items)) {
        return 0;
    }
    if (item_count == 0) {
        success = knapsack_result_replace(result, 0, NULL, 0);
        goto cleanup;
    }
    if (item_count == SIZE_MAX ||
        item_count + 1 > SIZE_MAX / sizeof(*suffix_values) ||
        item_count + 1 > SIZE_MAX / sizeof(*maximum_density_item)) {
        goto cleanup;
    }
    suffix_values =
        (int64_t *)calloc(item_count + 1, sizeof(*suffix_values));
    suffix_zero_weight_values = (int64_t *)calloc(
        item_count + 1, sizeof(*suffix_zero_weight_values));
    maximum_density_item = (size_t *)malloc(
        (item_count + 1) * sizeof(*maximum_density_item));
    if (suffix_values == NULL || suffix_zero_weight_values == NULL ||
        maximum_density_item == NULL) {
        goto cleanup;
    }
    maximum_density_item[item_count] = SIZE_MAX;
    for (i = item_count; i > 0; --i) {
        size_t item = i - 1;

        suffix_values[item] =
            add_bound_saturated(suffix_values[item + 1],
                                (uint64_t)items[item].value);
        suffix_zero_weight_values[item] =
            suffix_zero_weight_values[item + 1];
        if (items[item].weight == 0) {
            suffix_zero_weight_values[item] =
                add_bound_saturated(suffix_zero_weight_values[item],
                                    (uint64_t)items[item].value);
        }
        maximum_density_item[item] = maximum_density_item[item + 1];
        if (items[item].weight != 0 && items[item].value > 0) {
            maximum_density_item[item] = item;
        }
    }

    root.level = 0;
    root.weight = 0;
    root.value = 0;
    root.bound = density_upper_bound(
        items, suffix_values, suffix_zero_weight_values,
        maximum_density_item, 0, 0, 0, capacity);
    root.parent = SIZE_MAX;
    root.original_index = SIZE_MAX;
    root.included = 0;
    if (!bound_node_append(&nodes, &node_count, &node_capacity, &root,
                           &best_node)) {
        goto cleanup;
    }
    if (root.bound > best_value &&
        !open_set_push(&open_set, nodes, best_node)) {
        goto cleanup;
    }

    while (!open_set_is_empty(&open_set)) {
        size_t node_index = open_set_pop(&open_set, nodes);
        BoundNode node = nodes[node_index];

        if (node.value > best_value) {
            best_value = node.value;
            best_node = node_index;
        }
        if (bound_cannot_improve(node.bound, best_value) ||
            node.level == item_count) {
            continue;
        }

        if (policy == OPEN_LIFO) {
            if (!branch_bound_add_child(
                    items, suffix_values,
                    suffix_zero_weight_values, maximum_density_item,
                    capacity, &node, node_index,
                    node.level, 0, best_value, &nodes, &node_count,
                    &node_capacity, &open_set) ||
                !branch_bound_add_child(
                    items, suffix_values,
                    suffix_zero_weight_values, maximum_density_item,
                    capacity, &node, node_index,
                    node.level, 1, best_value, &nodes, &node_count,
                    &node_capacity, &open_set)) {
                goto cleanup;
            }
        } else {
            if (!branch_bound_add_child(
                    items, suffix_values,
                    suffix_zero_weight_values, maximum_density_item,
                    capacity, &node, node_index,
                    node.level, 1, best_value, &nodes, &node_count,
                    &node_capacity, &open_set) ||
                !branch_bound_add_child(
                    items, suffix_values,
                    suffix_zero_weight_values, maximum_density_item,
                    capacity, &node, node_index,
                    node.level, 0, best_value, &nodes, &node_count,
                    &node_capacity, &open_set)) {
                goto cleanup;
            }
        }
    }

    {
        size_t current = best_node;

        while (current != SIZE_MAX) {
            if (nodes[current].included) {
                if (selected_count == SIZE_MAX) {
                    goto cleanup;
                }
                ++selected_count;
            }
            current = nodes[current].parent;
        }
    }
    if (selected_count != 0) {
        size_t write = selected_count;
        size_t current = best_node;

        if (selected_count > SIZE_MAX / sizeof(*selected_indices)) {
            goto cleanup;
        }
        selected_indices = (size_t *)malloc(
            selected_count * sizeof(*selected_indices));
        if (selected_indices == NULL) {
            goto cleanup;
        }
        while (current != SIZE_MAX) {
            if (nodes[current].included) {
                selected_indices[--write] = nodes[current].original_index;
            }
            current = nodes[current].parent;
        }
    }
    success = knapsack_result_replace(result, best_value, selected_indices,
                                      selected_count);

cleanup:
    free(selected_indices);
    free(open_set.indices);
    free(nodes);
    free(maximum_density_item);
    free(suffix_zero_weight_values);
    free(suffix_values);
    free(items);
    return success;
}

int branch_bound_fifo_knapsack_01(const size_t *weights,
                                  const int64_t *values,
                                  size_t item_count,
                                  size_t capacity,
                                  KnapsackSearchResult *result)
{
    return branch_bound_knapsack(weights, values, item_count, capacity,
                                 OPEN_FIFO, result);
}

int branch_bound_lifo_knapsack_01(const size_t *weights,
                                  const int64_t *values,
                                  size_t item_count,
                                  size_t capacity,
                                  KnapsackSearchResult *result)
{
    return branch_bound_knapsack(weights, values, item_count, capacity,
                                 OPEN_LIFO, result);
}

int branch_bound_least_cost_knapsack_01(const size_t *weights,
                                        const int64_t *values,
                                        size_t item_count,
                                        size_t capacity,
                                        KnapsackSearchResult *result)
{
    return branch_bound_knapsack(weights, values, item_count, capacity,
                                 OPEN_BEST_BOUND, result);
}
