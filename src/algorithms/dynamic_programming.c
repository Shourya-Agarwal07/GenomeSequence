#include "genome/algorithms/dynamic_programming.h"

#include <limits.h>
#include <stdlib.h>
#include <string.h>

static int knapsack_01_compute(const size_t *weights,
                               const int64_t *values,
                               size_t item_count,
                               size_t capacity,
                               int64_t *maximum_value,
                               size_t **selected_indices,
                               size_t *selected_count)
{
    int64_t *table = NULL;
    unsigned char *take = NULL;
    size_t *indices = NULL;
    size_t selection_size = 0;
    size_t write_index;
    size_t columns;
    size_t i;
    int reconstruct = selected_indices != NULL;
    int success = 0;

    if (maximum_value != NULL) {
        *maximum_value = 0;
    }
    if (selected_indices != NULL) {
        *selected_indices = NULL;
    }
    if (selected_count != NULL) {
        *selected_count = 0;
    }
    if (maximum_value == NULL ||
        (reconstruct != (selected_count != NULL))) {
        return 0;
    }
    if (item_count == 0) {
        return 1;
    }
    if (weights == NULL || values == NULL) {
        return 0;
    }

    for (i = 0; i < item_count; ++i) {
        if (values[i] < 0) {
            return 0;
        }
    }
    if (capacity > SIZE_MAX / sizeof(*table) - 1) {
        return 0;
    }

    columns = capacity + 1;
    if (reconstruct && item_count > SIZE_MAX / columns) {
        return 0;
    }

    table = (int64_t *)calloc(capacity + 1, sizeof(*table));
    if (table == NULL) {
        return 0;
    }
    if (reconstruct) {
        size_t decision_cells = item_count * columns;

        if (decision_cells > SIZE_MAX / sizeof(*take)) {
            goto cleanup;
        }
        take = (unsigned char *)calloc(decision_cells, sizeof(*take));
        if (take == NULL) {
            goto cleanup;
        }
    }

    for (i = 0; i < item_count; ++i) {
        size_t current_capacity;

        if (weights[i] > capacity) {
            continue;
        }

        current_capacity = capacity;
        for (;;) {
            size_t previous_capacity = current_capacity - weights[i];
            int64_t candidate;

            if (values[i] > INT64_MAX - table[previous_capacity]) {
                goto cleanup;
            }
            candidate = table[previous_capacity] + values[i];
            if (candidate > table[current_capacity]) {
                table[current_capacity] = candidate;
                if (reconstruct) {
                    take[i * columns + current_capacity] = 1;
                }
            }

            if (current_capacity == weights[i]) {
                break;
            }
            --current_capacity;
        }
    }

    if (reconstruct) {
        size_t remaining_capacity = capacity;

        for (i = item_count; i > 0; --i) {
            size_t item = i - 1;

            if (take[item * columns + remaining_capacity]) {
                if (weights[item] > remaining_capacity) {
                    goto cleanup;
                }
                ++selection_size;
                remaining_capacity -= weights[item];
            }
        }

        if (selection_size > SIZE_MAX / sizeof(*indices)) {
            goto cleanup;
        }
        if (selection_size != 0) {
            indices = (size_t *)malloc(selection_size * sizeof(*indices));
            if (indices == NULL) {
                goto cleanup;
            }
        }

        remaining_capacity = capacity;
        write_index = selection_size;
        for (i = item_count; i > 0; --i) {
            size_t item = i - 1;

            if (take[item * columns + remaining_capacity]) {
                indices[--write_index] = item;
                remaining_capacity -= weights[item];
            }
        }
    }

    *maximum_value = table[capacity];
    if (reconstruct) {
        *selected_indices = indices;
        *selected_count = selection_size;
        indices = NULL;
    }
    success = 1;

cleanup:
    free(indices);
    free(take);
    free(table);
    return success;
}

int dp_knapsack_01(const size_t *weights,
                   const int64_t *values,
                   size_t item_count,
                   size_t capacity,
                   int64_t *maximum_value)
{
    return knapsack_01_compute(weights, values, item_count, capacity,
                               maximum_value, NULL, NULL);
}

int dp_knapsack_01_with_selection(const size_t *weights,
                                  const int64_t *values,
                                  size_t item_count,
                                  size_t capacity,
                                  int64_t *maximum_value,
                                  size_t **selected_indices,
                                  size_t *selected_count)
{
    return knapsack_01_compute(weights, values, item_count, capacity,
                               maximum_value, selected_indices,
                               selected_count);
}

int dp_longest_common_subsequence_n(const char *first,
                                    size_t first_length,
                                    const char *second,
                                    size_t second_length,
                                    size_t *lcs_length,
                                    char **lcs)
{
    size_t *table = NULL;
    char *result = NULL;
    size_t rows;
    size_t columns;
    size_t cell_count;
    size_t length;
    size_t i;
    size_t j;
    int success = 0;

    if (lcs_length != NULL) {
        *lcs_length = 0;
    }
    if (lcs != NULL) {
        *lcs = NULL;
    }
    if (first == NULL || second == NULL ||
        lcs_length == NULL || lcs == NULL ||
        first_length == SIZE_MAX || second_length == SIZE_MAX) {
        return 0;
    }

    rows = first_length + 1;
    columns = second_length + 1;
    if (rows > SIZE_MAX / columns) {
        return 0;
    }
    cell_count = rows * columns;
    if (cell_count > SIZE_MAX / sizeof(*table)) {
        return 0;
    }

    for (i = 0; i < first_length; ++i) {
        if (first[i] == '\0') {
            return 0;
        }
    }
    for (j = 0; j < second_length; ++j) {
        if (second[j] == '\0') {
            return 0;
        }
    }
    if (first_length == 0 || second_length == 0) {
        result = (char *)malloc(1);
        if (result == NULL) {
            return 0;
        }
        result[0] = '\0';
        *lcs = result;
        return 1;
    }

    table = (size_t *)calloc(cell_count, sizeof(*table));
    if (table == NULL) {
        return 0;
    }

    for (i = 1; i < rows; ++i) {
        for (j = 1; j < columns; ++j) {
            size_t current = i * columns + j;

            if (first[i - 1] == second[j - 1]) {
                table[current] = table[(i - 1) * columns + j - 1] + 1;
            } else {
                size_t above = table[(i - 1) * columns + j];
                size_t left = table[i * columns + j - 1];
                table[current] = above > left ? above : left;
            }
        }
    }

    length = table[first_length * columns + second_length];
    if (length == SIZE_MAX) {
        goto cleanup;
    }
    result = (char *)malloc(length + 1);
    if (result == NULL) {
        goto cleanup;
    }

    i = first_length;
    j = second_length;
    while (i != 0 && j != 0) {
        if (first[i - 1] == second[j - 1]) {
            result[--length] = first[i - 1];
            --i;
            --j;
        } else if (table[(i - 1) * columns + j] >=
                   table[i * columns + j - 1]) {
            --i;
        } else {
            --j;
        }
    }

    length = table[first_length * columns + second_length];
    result[length] = '\0';
    *lcs_length = length;
    *lcs = result;
    result = NULL;
    success = 1;

cleanup:
    free(result);
    free(table);
    return success;
}

int dp_longest_common_subsequence(const char *first,
                                  const char *second,
                                  size_t *lcs_length,
                                  char **lcs)
{
    if (lcs_length != NULL) {
        *lcs_length = 0;
    }
    if (lcs != NULL) {
        *lcs = NULL;
    }
    if (first == NULL || second == NULL ||
        lcs_length == NULL || lcs == NULL) {
        return 0;
    }

    return dp_longest_common_subsequence_n(
        first, strlen(first), second, strlen(second), lcs_length, lcs);
}

static int multiply_u64_checked(uint64_t left,
                                uint64_t right,
                                uint64_t *product)
{
    if (left != 0 && right > UINT64_MAX / left) {
        return 0;
    }
    *product = left * right;
    return 1;
}

static size_t decimal_digit_count(size_t value)
{
    size_t digits = 1;

    while (value >= 10) {
        value /= 10;
        ++digits;
    }
    return digits;
}

static int matrix_chain_output_size(size_t matrix_count,
                                   size_t *allocation_size)
{
    size_t string_length = 0;
    size_t internal_count = matrix_count - 1;
    size_t overhead;
    size_t matrix;

    for (matrix = 1; matrix <= matrix_count; ++matrix) {
        size_t label_length = 1 + decimal_digit_count(matrix);

        if (label_length > SIZE_MAX - string_length) {
            return 0;
        }
        string_length += label_length;
    }

    if (internal_count > SIZE_MAX / 5) {
        return 0;
    }
    overhead = internal_count * 5;
    if (overhead > SIZE_MAX - string_length) {
        return 0;
    }
    string_length += overhead;
    if (string_length == SIZE_MAX) {
        return 0;
    }
    *allocation_size = string_length + 1;
    return 1;
}

static void write_matrix_label(size_t matrix_number,
                               char *output,
                               size_t *position)
{
    char digits[sizeof(size_t) * CHAR_BIT];
    size_t digit_count = 0;

    output[(*position)++] = 'A';
    do {
        digits[digit_count++] = (char)('0' + matrix_number % 10);
        matrix_number /= 10;
    } while (matrix_number != 0);

    while (digit_count != 0) {
        output[(*position)++] = digits[--digit_count];
    }
}

static void write_matrix_chain_order(size_t first,
                                     size_t last,
                                     size_t matrix_count,
                                     const size_t *split_table,
                                     char *output,
                                     size_t *position)
{
    if (first == last) {
        write_matrix_label(first + 1, output, position);
        return;
    }

    {
        size_t split = split_table[first * matrix_count + last];

        output[(*position)++] = '(';
        write_matrix_chain_order(first, split, matrix_count, split_table,
                                 output, position);
        output[(*position)++] = ' ';
        output[(*position)++] = 'x';
        output[(*position)++] = ' ';
        write_matrix_chain_order(split + 1, last, matrix_count, split_table,
                                 output, position);
        output[(*position)++] = ')';
    }
}

int dp_matrix_chain_multiplication(const size_t *dimensions,
                                   size_t dimension_count,
                                   uint64_t *minimum_cost,
                                   char **parenthesization)
{
    uint64_t *cost_table = NULL;
    unsigned char *cost_valid = NULL;
    size_t *split_table = NULL;
    char *order = NULL;
    size_t matrix_count;
    size_t cell_count;
    size_t order_allocation_size;
    size_t i;
    size_t chain_length;
    uint64_t result_cost;
    size_t output_position = 0;
    int success = 0;

    if (minimum_cost != NULL) {
        *minimum_cost = 0;
    }
    if (parenthesization != NULL) {
        *parenthesization = NULL;
    }
    if (dimensions == NULL || minimum_cost == NULL ||
        parenthesization == NULL || dimension_count < 2) {
        return 0;
    }

    matrix_count = dimension_count - 1;
    if (matrix_count > SIZE_MAX / matrix_count) {
        return 0;
    }
    cell_count = matrix_count * matrix_count;
    if (cell_count > SIZE_MAX / sizeof(*cost_table) ||
        cell_count > SIZE_MAX / sizeof(*split_table) ||
        cell_count > SIZE_MAX / sizeof(*cost_valid) ||
        !matrix_chain_output_size(matrix_count, &order_allocation_size)) {
        return 0;
    }

    for (i = 0; i < dimension_count; ++i) {
        uint64_t converted_dimension = (uint64_t)dimensions[i];

        if (dimensions[i] == 0 ||
            (size_t)converted_dimension != dimensions[i]) {
            return 0;
        }
    }

    cost_table = (uint64_t *)calloc(cell_count, sizeof(*cost_table));
    if (cost_table == NULL) {
        goto cleanup;
    }
    cost_valid = (unsigned char *)calloc(cell_count, sizeof(*cost_valid));
    if (cost_valid == NULL) {
        goto cleanup;
    }
    split_table = (size_t *)malloc(cell_count * sizeof(*split_table));
    if (split_table == NULL) {
        goto cleanup;
    }
    for (i = 0; i < matrix_count; ++i) {
        cost_valid[i * matrix_count + i] = 1;
    }

    for (chain_length = 2; chain_length <= matrix_count; ++chain_length) {
        size_t first;

        for (first = 0; first <= matrix_count - chain_length; ++first) {
            size_t last = first + chain_length - 1;
            size_t split;
            size_t cell = first * matrix_count + last;
            int found_cost = 0;

            for (split = first; split < last; ++split) {
                uint64_t left_product;
                uint64_t multiplication_cost;
                uint64_t subchain_cost;
                uint64_t candidate_cost;
                size_t left_cell = first * matrix_count + split;
                size_t right_cell = (split + 1) * matrix_count + last;
                uint64_t left_cost = cost_table[left_cell];
                uint64_t right_cost = cost_table[right_cell];

                if (!cost_valid[left_cell] || !cost_valid[right_cell] ||
                    !multiply_u64_checked((uint64_t)dimensions[first],
                                          (uint64_t)dimensions[split + 1],
                                          &left_product) ||
                    !multiply_u64_checked(
                        left_product, (uint64_t)dimensions[last + 1],
                        &multiplication_cost) ||
                    left_cost > UINT64_MAX - right_cost) {
                    continue;
                }
                subchain_cost = left_cost + right_cost;
                if (multiplication_cost > UINT64_MAX - subchain_cost) {
                    continue;
                }
                candidate_cost = subchain_cost + multiplication_cost;

                if (!found_cost || candidate_cost < cost_table[cell]) {
                    cost_table[cell] = candidate_cost;
                    split_table[cell] = split;
                    found_cost = 1;
                }
            }
            if (found_cost) {
                cost_valid[cell] = 1;
            }
        }
    }

    if (!cost_valid[matrix_count - 1]) {
        goto cleanup;
    }

    order = (char *)malloc(order_allocation_size);
    if (order == NULL) {
        goto cleanup;
    }
    write_matrix_chain_order(0, matrix_count - 1, matrix_count, split_table,
                             order, &output_position);
    if (output_position != order_allocation_size - 1) {
        goto cleanup;
    }
    order[output_position] = '\0';

    result_cost = cost_table[matrix_count - 1];
    *minimum_cost = result_cost;
    *parenthesization = order;
    order = NULL;
    success = 1;

cleanup:
    free(order);
    free(split_table);
    free(cost_valid);
    free(cost_table);
    return success;
}

static int add_u64_checked(uint64_t left, uint64_t right, uint64_t *sum)
{
    if (left > UINT64_MAX - right) {
        return 0;
    }
    *sum = left + right;
    return 1;
}

void dp_optimal_bst_destroy(DpOptimalBstNode *root)
{
    if (root == NULL) {
        return;
    }
    dp_optimal_bst_destroy(root->left);
    dp_optimal_bst_destroy(root->right);
    free(root);
}

static int build_optimal_bst_tree(size_t first,
                                  size_t last,
                                  size_t key_count,
                                  const int *keys,
                                  const uint64_t *frequencies,
                                  const size_t *root_table,
                                  DpOptimalBstNode **tree)
{
    DpOptimalBstNode *node;
    size_t root_index;

    *tree = NULL;
    if (first == last) {
        return 1;
    }

    root_index = root_table[first * key_count + last - 1];
    node = (DpOptimalBstNode *)malloc(sizeof(*node));
    if (node == NULL) {
        return 0;
    }
    node->key = keys[root_index];
    node->frequency = frequencies[root_index];
    node->left = NULL;
    node->right = NULL;

    if (!build_optimal_bst_tree(first, root_index, key_count, keys,
                                frequencies, root_table, &node->left) ||
        !build_optimal_bst_tree(root_index + 1, last, key_count, keys,
                                frequencies, root_table, &node->right)) {
        dp_optimal_bst_destroy(node);
        return 0;
    }

    *tree = node;
    return 1;
}

int dp_optimal_bst(const int *keys,
                   const uint64_t *frequencies,
                   size_t key_count,
                   uint64_t *minimum_cost,
                   DpOptimalBstNode **root)
{
    uint64_t *cost_table = NULL;
    size_t *root_table = NULL;
    unsigned char *cost_valid = NULL;
    uint64_t *frequency_prefix = NULL;
    DpOptimalBstNode *tree = NULL;
    size_t cell_count;
    size_t prefix_count;
    size_t i;
    size_t interval_length;
    size_t result_cell;
    int success = 0;

    if (minimum_cost != NULL) {
        *minimum_cost = 0;
    }
    if (root != NULL) {
        *root = NULL;
    }
    if (minimum_cost == NULL || root == NULL) {
        return 0;
    }
    if (key_count == 0) {
        return 1;
    }
    if (keys == NULL || frequencies == NULL || key_count == SIZE_MAX ||
        key_count > SIZE_MAX / key_count) {
        return 0;
    }

    cell_count = key_count * key_count;
    prefix_count = key_count + 1;
    if (cell_count > SIZE_MAX / sizeof(*cost_table) ||
        cell_count > SIZE_MAX / sizeof(*root_table) ||
        cell_count > SIZE_MAX / sizeof(*cost_valid) ||
        prefix_count > SIZE_MAX / sizeof(*frequency_prefix)) {
        return 0;
    }

    frequency_prefix =
        (uint64_t *)malloc(prefix_count * sizeof(*frequency_prefix));
    if (frequency_prefix == NULL) {
        goto cleanup;
    }
    frequency_prefix[0] = 0;
    for (i = 0; i < key_count; ++i) {
        if ((i != 0 && keys[i - 1] >= keys[i]) ||
            !add_u64_checked(frequency_prefix[i], frequencies[i],
                             &frequency_prefix[i + 1])) {
            goto cleanup;
        }
    }

    cost_table = (uint64_t *)calloc(cell_count, sizeof(*cost_table));
    if (cost_table == NULL) {
        goto cleanup;
    }
    root_table = (size_t *)malloc(cell_count * sizeof(*root_table));
    if (root_table == NULL) {
        goto cleanup;
    }
    cost_valid = (unsigned char *)calloc(cell_count, sizeof(*cost_valid));
    if (cost_valid == NULL) {
        goto cleanup;
    }

    for (interval_length = 1; interval_length <= key_count;
         ++interval_length) {
        size_t first;

        for (first = 0; first <= key_count - interval_length; ++first) {
            size_t last = first + interval_length;
            size_t cell = first * key_count + last - 1;
            size_t candidate_root;
            uint64_t interval_frequency =
                frequency_prefix[last] - frequency_prefix[first];
            int found_cost = 0;

            for (candidate_root = first; candidate_root < last;
                 ++candidate_root) {
                uint64_t left_cost = 0;
                uint64_t right_cost = 0;
                uint64_t child_cost;
                uint64_t candidate_cost;
                size_t left_cell = 0;
                size_t right_cell = 0;

                if (candidate_root != first) {
                    left_cell =
                        first * key_count + candidate_root - 1;
                    if (!cost_valid[left_cell]) {
                        continue;
                    }
                    left_cost = cost_table[left_cell];
                }
                if (candidate_root + 1 != last) {
                    right_cell =
                        (candidate_root + 1) * key_count + last - 1;
                    if (!cost_valid[right_cell]) {
                        continue;
                    }
                    right_cost = cost_table[right_cell];
                }
                if (!add_u64_checked(left_cost, right_cost, &child_cost) ||
                    !add_u64_checked(child_cost, interval_frequency,
                                     &candidate_cost)) {
                    continue;
                }

                if (!found_cost || candidate_cost < cost_table[cell]) {
                    cost_table[cell] = candidate_cost;
                    root_table[cell] = candidate_root;
                    found_cost = 1;
                }
            }
            if (found_cost) {
                cost_valid[cell] = 1;
            }
        }
    }

    result_cell = key_count - 1;
    if (!cost_valid[result_cell] ||
        !build_optimal_bst_tree(0, key_count, key_count, keys, frequencies,
                                root_table, &tree)) {
        goto cleanup;
    }

    *minimum_cost = cost_table[result_cell];
    *root = tree;
    tree = NULL;
    success = 1;

cleanup:
    dp_optimal_bst_destroy(tree);
    free(cost_valid);
    free(root_table);
    free(cost_table);
    free(frequency_prefix);
    return success;
}

static int add_i64_checked(int64_t left, int64_t right, int64_t *sum)
{
    if ((right > 0 && left > INT64_MAX - right) ||
        (right < 0 && left < INT64_MIN - right)) {
        return 0;
    }
    *sum = left + right;
    return 1;
}

void dp_bellman_ford_result_destroy(DpBellmanFordResult *result)
{
    if (result == NULL) {
        return;
    }
    free(result->distances);
    free(result->predecessors);
    free(result->reachable);
    result->vertex_count = 0;
    result->source = 0;
    result->distances = NULL;
    result->predecessors = NULL;
    result->reachable = NULL;
    result->negative_cycle_detected = 0;
}

int dp_bellman_ford(size_t vertex_count,
                    const DpBellmanFordEdge *edges,
                    size_t edge_count,
                    size_t source,
                    DpBellmanFordResult *result)
{
    int64_t *distances = NULL;
    size_t *predecessors = NULL;
    unsigned char *reachable = NULL;
    size_t i;
    size_t pass;
    int negative_cycle_detected = 0;
    int success = 0;

    if (result == NULL || vertex_count == 0 || source >= vertex_count ||
        (edge_count != 0 && edges == NULL) ||
        edge_count > SIZE_MAX / sizeof(*edges) ||
        vertex_count > SIZE_MAX / sizeof(*distances) ||
        vertex_count > SIZE_MAX / sizeof(*predecessors) ||
        vertex_count > SIZE_MAX / sizeof(*reachable)) {
        return 0;
    }
    for (i = 0; i < edge_count; ++i) {
        if (edges[i].from >= vertex_count || edges[i].to >= vertex_count) {
            return 0;
        }
    }

    distances = (int64_t *)calloc(vertex_count, sizeof(*distances));
    predecessors = (size_t *)malloc(vertex_count * sizeof(*predecessors));
    reachable = (unsigned char *)calloc(vertex_count, sizeof(*reachable));
    if (distances == NULL || predecessors == NULL || reachable == NULL) {
        goto cleanup;
    }
    for (i = 0; i < vertex_count; ++i) {
        predecessors[i] = SIZE_MAX;
    }
    reachable[source] = 1;

    for (pass = 0; pass < vertex_count - 1; ++pass) {
        int changed = 0;
        size_t edge_index;

        for (edge_index = 0; edge_index < edge_count; ++edge_index) {
            const DpBellmanFordEdge *edge = &edges[edge_index];
            int64_t candidate;

            if (!reachable[edge->from]) {
                continue;
            }
            if (!add_i64_checked(distances[edge->from], edge->weight,
                                 &candidate)) {
                goto cleanup;
            }
            if (!reachable[edge->to] || candidate < distances[edge->to]) {
                distances[edge->to] = candidate;
                predecessors[edge->to] = edge->from;
                reachable[edge->to] = 1;
                changed = 1;
            }
        }
        if (!changed) {
            break;
        }
    }

    for (i = 0; i < edge_count; ++i) {
        const DpBellmanFordEdge *edge = &edges[i];
        int64_t candidate;

        if (!reachable[edge->from]) {
            continue;
        }
        if (!add_i64_checked(distances[edge->from], edge->weight,
                             &candidate)) {
            goto cleanup;
        }
        if (!reachable[edge->to] || candidate < distances[edge->to]) {
            negative_cycle_detected = 1;
            break;
        }
    }

    dp_bellman_ford_result_destroy(result);
    result->vertex_count = vertex_count;
    result->source = source;
    result->distances = distances;
    result->predecessors = predecessors;
    result->reachable = reachable;
    result->negative_cycle_detected = negative_cycle_detected;
    distances = NULL;
    predecessors = NULL;
    reachable = NULL;
    success = 1;

cleanup:
    free(reachable);
    free(predecessors);
    free(distances);
    return success;
}

int dp_bellman_ford_reconstruct_path(const DpBellmanFordResult *result,
                                     size_t destination,
                                     size_t **path,
                                     size_t *path_length)
{
    size_t *reconstructed = NULL;
    size_t reconstructed_length = 0;
    size_t current;
    size_t i;
    int success = 0;

    if (path != NULL) {
        *path = NULL;
    }
    if (path_length != NULL) {
        *path_length = 0;
    }
    if (result == NULL || path == NULL || path_length == NULL ||
        result->vertex_count == 0 || result->source >= result->vertex_count ||
        destination >= result->vertex_count || result->distances == NULL ||
        result->predecessors == NULL || result->reachable == NULL ||
        result->negative_cycle_detected || !result->reachable[destination] ||
        result->vertex_count > SIZE_MAX / sizeof(*reconstructed)) {
        return 0;
    }

    reconstructed =
        (size_t *)malloc(result->vertex_count * sizeof(*reconstructed));
    if (reconstructed == NULL) {
        return 0;
    }

    current = destination;
    while (reconstructed_length < result->vertex_count) {
        if (current >= result->vertex_count) {
            goto cleanup;
        }
        reconstructed[reconstructed_length++] = current;
        if (current == result->source) {
            break;
        }
        current = result->predecessors[current];
        if (current == SIZE_MAX) {
            goto cleanup;
        }
    }
    if (reconstructed[reconstructed_length - 1] != result->source) {
        goto cleanup;
    }

    for (i = 0; i < reconstructed_length / 2; ++i) {
        size_t opposite = reconstructed_length - 1 - i;
        size_t temporary = reconstructed[i];

        reconstructed[i] = reconstructed[opposite];
        reconstructed[opposite] = temporary;
    }

    *path = reconstructed;
    *path_length = reconstructed_length;
    reconstructed = NULL;
    success = 1;

cleanup:
    free(reconstructed);
    return success;
}

static int weighted_matrix_cell_count(const DpWeightedDirectedMatrix *graph,
                                      size_t *cell_count)
{
    size_t vertices;

    if (graph == NULL || cell_count == NULL || graph->vertex_count == 0 ||
        graph->weights == NULL || graph->edge_present == NULL) {
        return 0;
    }
    vertices = graph->vertex_count;
    if (vertices > SIZE_MAX / vertices) {
        return 0;
    }
    *cell_count = vertices * vertices;
    return *cell_count <= SIZE_MAX / sizeof(*graph->weights);
}

void dp_floyd_warshall_result_destroy(DpFloydWarshallResult *result)
{
    if (result == NULL) {
        return;
    }
    free(result->distances);
    free(result->next_hops);
    free(result->reachable);
    free(result->negative_cycle_affected);
    result->vertex_count = 0;
    result->distances = NULL;
    result->next_hops = NULL;
    result->reachable = NULL;
    result->negative_cycle_affected = NULL;
    result->negative_cycle_detected = 0;
}

int dp_floyd_warshall(const DpWeightedDirectedMatrix *graph,
                      DpFloydWarshallResult *result)
{
    int64_t *distances = NULL;
    size_t *next_hops = NULL;
    unsigned char *reachable = NULL;
    unsigned char *negative_cycle_affected = NULL;
    size_t cell_count;
    size_t vertices;
    size_t i;
    size_t j;
    size_t k;
    int negative_cycle_detected = 0;
    int success = 0;

    if (result == NULL || !weighted_matrix_cell_count(graph, &cell_count)) {
        return 0;
    }
    vertices = graph->vertex_count;
    if (cell_count > SIZE_MAX / sizeof(*distances) ||
        cell_count > SIZE_MAX / sizeof(*next_hops) ||
        cell_count > SIZE_MAX / sizeof(*reachable) ||
        cell_count > SIZE_MAX / sizeof(*negative_cycle_affected)) {
        return 0;
    }

    distances = (int64_t *)malloc(cell_count * sizeof(*distances));
    next_hops = (size_t *)malloc(cell_count * sizeof(*next_hops));
    reachable =
        (unsigned char *)calloc(cell_count, sizeof(*reachable));
    negative_cycle_affected =
        (unsigned char *)calloc(cell_count, sizeof(*negative_cycle_affected));
    if (distances == NULL || next_hops == NULL || reachable == NULL ||
        negative_cycle_affected == NULL) {
        goto cleanup;
    }

    for (i = 0; i < vertices; ++i) {
        for (j = 0; j < vertices; ++j) {
            size_t cell = i * vertices + j;

            next_hops[cell] = SIZE_MAX;
            if (i == j) {
                distances[cell] = 0;
                reachable[cell] = 1;
                next_hops[cell] = i;
                if (graph->edge_present[cell] &&
                    graph->weights[cell] < 0) {
                    distances[cell] = graph->weights[cell];
                }
            } else if (graph->edge_present[cell]) {
                distances[cell] = graph->weights[cell];
                reachable[cell] = 1;
                next_hops[cell] = j;
            } else {
                distances[cell] = 0;
            }
        }
    }

    for (k = 0; k < vertices; ++k) {
        for (i = 0; i < vertices; ++i) {
            size_t to_intermediate = i * vertices + k;

            if (!reachable[to_intermediate]) {
                continue;
            }
            for (j = 0; j < vertices; ++j) {
                size_t from_intermediate = k * vertices + j;
                size_t cell = i * vertices + j;
                int64_t candidate;

                if (!reachable[from_intermediate]) {
                    continue;
                }
                if (!add_i64_checked(distances[to_intermediate],
                                     distances[from_intermediate],
                                     &candidate)) {
                    goto cleanup;
                }
                if (!reachable[cell] || candidate < distances[cell]) {
                    distances[cell] = candidate;
                    reachable[cell] = 1;
                    next_hops[cell] = next_hops[to_intermediate];
                }
            }
        }
    }

    for (k = 0; k < vertices; ++k) {
        if (distances[k * vertices + k] < 0) {
            negative_cycle_detected = 1;
            for (i = 0; i < vertices; ++i) {
                if (!reachable[i * vertices + k]) {
                    continue;
                }
                for (j = 0; j < vertices; ++j) {
                    size_t cell = i * vertices + j;

                    if (reachable[k * vertices + j]) {
                        negative_cycle_affected[cell] = 1;
                        next_hops[cell] = SIZE_MAX;
                    }
                }
            }
        }
    }

    dp_floyd_warshall_result_destroy(result);
    result->vertex_count = vertices;
    result->distances = distances;
    result->next_hops = next_hops;
    result->reachable = reachable;
    result->negative_cycle_affected = negative_cycle_affected;
    result->negative_cycle_detected = negative_cycle_detected;
    distances = NULL;
    next_hops = NULL;
    reachable = NULL;
    negative_cycle_affected = NULL;
    success = 1;

cleanup:
    free(negative_cycle_affected);
    free(reachable);
    free(next_hops);
    free(distances);
    return success;
}

int dp_floyd_warshall_reconstruct_path(
    const DpFloydWarshallResult *result,
    size_t source,
    size_t destination,
    size_t **path,
    size_t *path_length)
{
    size_t *reconstructed = NULL;
    size_t reconstructed_length = 0;
    size_t current;
    int success = 0;

    if (path != NULL) {
        *path = NULL;
    }
    if (path_length != NULL) {
        *path_length = 0;
    }
    if (result == NULL || path == NULL || path_length == NULL ||
        result->vertex_count == 0 || source >= result->vertex_count ||
        destination >= result->vertex_count || result->next_hops == NULL ||
        result->reachable == NULL ||
        result->negative_cycle_affected == NULL ||
        result->vertex_count > SIZE_MAX / result->vertex_count ||
        !result->reachable[source * result->vertex_count + destination] ||
        result->negative_cycle_affected[
            source * result->vertex_count + destination] ||
        result->vertex_count > SIZE_MAX / sizeof(*reconstructed)) {
        return 0;
    }

    reconstructed = (size_t *)malloc(
        result->vertex_count * sizeof(*reconstructed));
    if (reconstructed == NULL) {
        return 0;
    }

    current = source;
    reconstructed[reconstructed_length++] = source;
    while (current != destination &&
           reconstructed_length < result->vertex_count) {
        size_t next =
            result->next_hops[current * result->vertex_count + destination];

        if (next == SIZE_MAX || next >= result->vertex_count) {
            goto cleanup;
        }
        reconstructed[reconstructed_length++] = next;
        current = next;
    }
    if (current != destination) {
        goto cleanup;
    }

    *path = reconstructed;
    *path_length = reconstructed_length;
    reconstructed = NULL;
    success = 1;

cleanup:
    free(reconstructed);
    return success;
}

void dp_tsp_result_destroy(DpTspResult *result)
{
    if (result == NULL) {
        return;
    }
    free(result->tour);
    result->vertex_count = 0;
    result->minimum_cost = 0;
    result->tour = NULL;
    result->tour_length = 0;
}

DpTspStatus dp_tsp_bitmask(const DpWeightedDirectedMatrix *graph,
                           size_t start,
                           DpTspResult *result)
{
    int64_t *costs = NULL;
    uint8_t *parents = NULL;
    size_t *tour = NULL;
    size_t vertex_count;
    size_t matrix_cells;
    size_t state_count;
    size_t state_cells;
    size_t start_mask;
    size_t full_mask;
    size_t mask;
    size_t best_last = SIZE_MAX;
    size_t tour_length;
    size_t position;
    int64_t best_cost = 0;
    int found_tour = 0;
    DpTspStatus status = DP_TSP_ERROR;

    if (result == NULL ||
        !weighted_matrix_cell_count(graph, &matrix_cells)) {
        return DP_TSP_ERROR;
    }
    vertex_count = graph->vertex_count;
    if (vertex_count > DP_TSP_MAX_VERTICES || start >= vertex_count ||
        vertex_count >= sizeof(size_t) * CHAR_BIT) {
        return DP_TSP_ERROR;
    }
    state_count = (size_t)1 << vertex_count;
    if (state_count > SIZE_MAX / vertex_count) {
        return DP_TSP_ERROR;
    }
    state_cells = state_count * vertex_count;
    if (state_cells > SIZE_MAX / sizeof(*costs) ||
        state_cells > SIZE_MAX / sizeof(*parents)) {
        return DP_TSP_ERROR;
    }
    if (matrix_cells > SIZE_MAX / sizeof(*graph->weights) ||
        vertex_count == SIZE_MAX ||
        vertex_count + 1 > SIZE_MAX / sizeof(*tour)) {
        return DP_TSP_ERROR;
    }

    costs = (int64_t *)calloc(state_cells, sizeof(*costs));
    parents = (uint8_t *)malloc(state_cells * sizeof(*parents));
    if (costs == NULL || parents == NULL) {
        goto cleanup;
    }
    memset(parents, UINT8_MAX, state_cells * sizeof(*parents));

    start_mask = (size_t)1 << start;
    full_mask = state_count - 1;
    costs[start_mask * vertex_count + start] = 0;

    for (mask = 0; mask < state_count; ++mask) {
        size_t last;

        if ((mask & start_mask) == 0) {
            continue;
        }
        for (last = 0; last < vertex_count; ++last) {
            size_t state;
            int base_state;
            size_t next;

            if ((mask & ((size_t)1 << last)) == 0) {
                continue;
            }
            state = mask * vertex_count + last;
            base_state = mask == start_mask && last == start;
            if (!base_state && parents[state] == UINT8_MAX) {
                continue;
            }

            for (next = 0; next < vertex_count; ++next) {
                size_t next_bit = (size_t)1 << next;
                size_t next_mask;
                size_t next_state;
                size_t edge_cell;
                int64_t candidate;

                if ((mask & next_bit) != 0) {
                    continue;
                }
                edge_cell = last * vertex_count + next;
                if (!graph->edge_present[edge_cell]) {
                    continue;
                }
                if (!add_i64_checked(costs[state],
                                     graph->weights[edge_cell],
                                     &candidate)) {
                    goto cleanup;
                }
                next_mask = mask | next_bit;
                next_state = next_mask * vertex_count + next;
                if (parents[next_state] == UINT8_MAX ||
                    candidate < costs[next_state]) {
                    costs[next_state] = candidate;
                    parents[next_state] = (uint8_t)last;
                }
            }
        }
    }

    {
        size_t last;

        for (last = 0; last < vertex_count; ++last) {
            size_t state = full_mask * vertex_count + last;
            size_t edge_cell = last * vertex_count + start;
            int base_state = full_mask == start_mask && last == start;
            int64_t candidate;

            if (!base_state && parents[state] == UINT8_MAX) {
                continue;
            }
            if (!graph->edge_present[edge_cell]) {
                continue;
            }
            if (!add_i64_checked(costs[state], graph->weights[edge_cell],
                                 &candidate)) {
                goto cleanup;
            }
            if (!found_tour || candidate < best_cost) {
                best_cost = candidate;
                best_last = last;
                found_tour = 1;
            }
        }
    }

    if (!found_tour) {
        dp_tsp_result_destroy(result);
        status = DP_TSP_NO_TOUR;
        goto cleanup;
    }

    {
        size_t selected_last = best_last;

        if (vertex_count + 1 > SIZE_MAX / sizeof(*tour)) {
            goto cleanup;
        }
        tour_length = vertex_count + 1;
        tour = (size_t *)malloc(tour_length * sizeof(*tour));
        if (tour == NULL) {
            goto cleanup;
        }
        tour[vertex_count] = start;
        position = vertex_count;
        mask = full_mask;
        while (position > 1) {
            size_t state = mask * vertex_count + selected_last;
            uint8_t previous = parents[state];

            if (previous == UINT8_MAX || previous >= vertex_count) {
                goto cleanup;
            }
            tour[position - 1] = selected_last;
            mask ^= (size_t)1 << selected_last;
            selected_last = previous;
            --position;
        }
        if (selected_last != start || mask != start_mask) {
            goto cleanup;
        }
        tour[0] = start;
    }

    dp_tsp_result_destroy(result);
    result->vertex_count = vertex_count;
    result->minimum_cost = best_cost;
    result->tour = tour;
    result->tour_length = tour_length;
    tour = NULL;
    status = DP_TSP_TOUR_FOUND;

cleanup:
    free(tour);
    free(parents);
    free(costs);
    return status;
}
