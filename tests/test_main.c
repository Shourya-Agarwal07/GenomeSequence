#include "genome/domain/dna_sequence.h"
#include "genome/io/fasta_reader.h"
#include "genome/algorithms/sequence_algorithms.h"
#include "genome/algorithms/divide_conquer.h"
#include "genome/algorithms/dynamic_programming.h"
#include "genome/algorithms/dna_alignment.h"
#include "genome/algorithms/backtracking_algorithms.h"
#include "genome/algorithms/greedy_algorithms.h"
#include "genome/graph/ford_fulkerson.h"
#include "genome/graph/bipartite_matching.h"
#include "genome/string_matching/string_matching.h"
#include "genome/benchmark/benchmark.h"
#include "genome/export/result_export.h"
#include "genome/sequence/sequence_statistics.h"
#include "genome/cli/cli.h"

#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef int (*TestFunction)(void);
typedef int (*SortFunction)(DnaSequence *,
                            size_t,
                            SequenceSortCriterion,
                            SequenceAlgorithmStats *);
typedef SequenceSearchResult (*SearchFunction)(const DnaSequence *,
                                              size_t,
                                              const DnaSequence *,
                                              SequenceSortCriterion,
                                              size_t *,
                                              SequenceAlgorithmStats *);

static int matrices_are_close(const DcaMatrix *left, const DcaMatrix *right)
{
    size_t count;
    size_t i;

    if (left->rows != right->rows || left->columns != right->columns) {
        return 0;
    }
    count = left->rows * left->columns;
    for (i = 0; i < count; ++i) {
        double difference = left->values[i] - right->values[i];
        double tolerance = 1e-10 * (1.0 + fabs(right->values[i]));

        if (!isfinite(left->values[i]) || !isfinite(right->values[i]) ||
            fabs(difference) > tolerance) {
            fprintf(stderr,
                    "Matrix mismatch at (%lu, %lu): %.17g != %.17g\n",
                    (unsigned long)(i / left->columns),
                    (unsigned long)(i % left->columns),
                    left->values[i], right->values[i]);
            return 0;
        }
    }
    return 1;
}

static int matrix_matches_values(const DcaMatrix *matrix,
                                 size_t rows,
                                 size_t columns,
                                 const double *expected)
{
    size_t count;
    size_t i;

    if (matrix->rows != rows || matrix->columns != columns) {
        return 0;
    }
    count = rows * columns;
    for (i = 0; i < count; ++i) {
        if (matrix->values[i] != expected[i]) {
            return 0;
        }
    }
    return 1;
}

static int test_maximum_subarray(void)
{
    const int64_t normal[] = {-2, 1, -3, 4, -1, 2, 1, -5, 4};
    const int64_t negative[] = {-8, -2, -5};
    const int64_t zeroes[] = {0, 0, 0};
    const int64_t single[] = {-42};
    const int64_t overflow[] = {INT64_MAX, 1};
    MaximumSubarrayResult result;
    DcaOperationCounts operations;

    if (!dca_maximum_subarray(normal, 9, &result, &operations) ||
        result.sum != 6 || result.start_index != 3 || result.end_index != 6 ||
        operations.additions == 0 || operations.comparisons == 0) {
        return 0;
    }
    if (!dca_maximum_subarray(negative, 3, &result, NULL) ||
        result.sum != -2 || result.start_index != 1 || result.end_index != 1) {
        return 0;
    }
    if (!dca_maximum_subarray(zeroes, 3, &result, NULL) ||
        result.sum != 0 || result.start_index > result.end_index ||
        result.end_index >= 3) {
        return 0;
    }
    if (!dca_maximum_subarray(single, 1, &result, NULL) ||
        result.sum != -42 || result.start_index != 0 || result.end_index != 0) {
        return 0;
    }
    if (dca_maximum_subarray(overflow, 2, &result, NULL) ||
        dca_maximum_subarray(NULL, 1, &result, NULL) ||
        dca_maximum_subarray(single, 0, &result, NULL) ||
        dca_maximum_subarray(single, 1, NULL, NULL)) {
        return 0;
    }

    return 1;
}

static int test_dp_knapsack_01(void)
{
    const size_t weights[] = {10, 20, 30};
    const int64_t values[] = {60, 100, 120};
    const size_t zero_weight[] = {0, 2, 3};
    const int64_t zero_weight_values[] = {7, 8, 11};
    const size_t single_weight[] = {2};
    const int64_t single_value[] = {5};
    const size_t overflow_weights[] = {1, 1};
    const int64_t overflow_values[] = {INT64_MAX, 1};
    const int64_t negative_value[] = {-1};
    int64_t maximum_value = -1;

    if (!dp_knapsack_01(weights, values, 3, 50, &maximum_value) ||
        maximum_value != 220 ||
        !dp_knapsack_01(zero_weight, zero_weight_values, 3, 3,
                        &maximum_value) ||
        maximum_value != 18 ||
        !dp_knapsack_01(single_weight, single_value, 1, 6,
                        &maximum_value) ||
        maximum_value != 5 ||
        !dp_knapsack_01(single_weight, single_value, 1, 0,
                        &maximum_value) ||
        maximum_value != 0 ||
        !dp_knapsack_01(NULL, NULL, 0, SIZE_MAX, &maximum_value) ||
        maximum_value != 0) {
        return 0;
    }

    maximum_value = -1;
    if (dp_knapsack_01(NULL, values, 3, 50, &maximum_value) ||
        maximum_value != 0 ||
        dp_knapsack_01(single_weight, negative_value, 1, 2,
                       &maximum_value) ||
        maximum_value != 0 ||
        dp_knapsack_01(overflow_weights, overflow_values, 2, 2,
                       &maximum_value) ||
        maximum_value != 0 ||
        dp_knapsack_01(single_weight, single_value, 1, SIZE_MAX,
                       &maximum_value) ||
        maximum_value != 0 ||
        dp_knapsack_01(single_weight, single_value, 1, 2, NULL)) {
        return 0;
    }

    return 1;
}

static int knapsack_selection_matches(const size_t *weights,
                                      const int64_t *values,
                                      size_t item_count,
                                      size_t capacity,
                                      const size_t *selected_indices,
                                      size_t selected_count,
                                      int64_t maximum_value,
                                      size_t *selected_weight)
{
    size_t total_weight = 0;
    int64_t total_value = 0;
    size_t i;

    if (selected_count > item_count ||
        ((selected_count == 0) != (selected_indices == NULL))) {
        return 0;
    }

    for (i = 0; i < selected_count; ++i) {
        size_t item = selected_indices[i];
        size_t j;

        if (item >= item_count ||
            (i != 0 && selected_indices[i - 1] >= item) ||
            weights[item] > capacity - total_weight ||
            values[item] > INT64_MAX - total_value) {
            return 0;
        }
        for (j = 0; j < i; ++j) {
            if (selected_indices[j] == item) {
                return 0;
            }
        }
        total_weight += weights[item];
        total_value += values[item];
    }

    if (selected_weight != NULL) {
        *selected_weight = total_weight;
    }
    return total_value == maximum_value;
}

static int test_dp_knapsack_01_reconstruction(void)
{
    const size_t weights[] = {10, 20, 30};
    const int64_t values[] = {60, 100, 120};
    const size_t zero_weights[] = {0, 2, 3};
    const int64_t zero_values[] = {7, 8, 11};
    const size_t zero_capacity_weights[] = {0, 1};
    const int64_t zero_capacity_values[] = {9, 5};
    const size_t no_fit_weights[] = {3, 4};
    const int64_t no_fit_values[] = {10, 20};
    const size_t tied_weights[] = {1, 1};
    const int64_t tied_values[] = {5, 5};
    size_t *selected_indices = NULL;
    size_t selected_count = 0;
    size_t selected_weight = 0;
    int64_t maximum_value = -1;
    int valid = 0;

    if (!dp_knapsack_01_with_selection(
            weights, values, 3, 50, &maximum_value,
            &selected_indices, &selected_count) ||
        maximum_value != 220 || selected_count != 2 ||
        selected_indices == NULL || selected_indices[0] != 1 ||
        selected_indices[1] != 2 ||
        !knapsack_selection_matches(weights, values, 3, 50,
                                    selected_indices, selected_count,
                                    maximum_value, &selected_weight) ||
        selected_weight != 50) {
        goto cleanup;
    }
    free(selected_indices);
    selected_indices = NULL;

    if (!dp_knapsack_01_with_selection(
            zero_weights, zero_values, 3, 3, &maximum_value,
            &selected_indices, &selected_count) ||
        maximum_value != 18 || selected_count != 2 ||
        selected_indices == NULL || selected_indices[0] != 0 ||
        selected_indices[1] != 2 ||
        !knapsack_selection_matches(zero_weights, zero_values, 3, 3,
                                    selected_indices, selected_count,
                                    maximum_value, &selected_weight) ||
        selected_weight != 3) {
        goto cleanup;
    }
    free(selected_indices);
    selected_indices = NULL;

    if (!dp_knapsack_01_with_selection(
            NULL, NULL, 0, 20, &maximum_value,
            &selected_indices, &selected_count) ||
        maximum_value != 0 || selected_indices != NULL ||
        selected_count != 0 ||
        !knapsack_selection_matches(NULL, NULL, 0, 20, selected_indices,
                                    selected_count, maximum_value,
                                    &selected_weight) ||
        selected_weight != 0) {
        goto cleanup;
    }

    if (!dp_knapsack_01_with_selection(
            zero_capacity_weights, zero_capacity_values, 2, 0,
            &maximum_value, &selected_indices, &selected_count) ||
        maximum_value != 9 || selected_count != 1 ||
        selected_indices == NULL || selected_indices[0] != 0 ||
        !knapsack_selection_matches(zero_capacity_weights,
                                    zero_capacity_values, 2, 0,
                                    selected_indices, selected_count,
                                    maximum_value, &selected_weight) ||
        selected_weight != 0) {
        goto cleanup;
    }
    free(selected_indices);
    selected_indices = NULL;

    if (!dp_knapsack_01_with_selection(
            no_fit_weights, no_fit_values, 2, 2, &maximum_value,
            &selected_indices, &selected_count) ||
        maximum_value != 0 || selected_indices != NULL ||
        selected_count != 0 ||
        !knapsack_selection_matches(no_fit_weights, no_fit_values, 2, 2,
                                    selected_indices, selected_count,
                                    maximum_value, &selected_weight) ||
        selected_weight != 0) {
        goto cleanup;
    }

    if (!dp_knapsack_01_with_selection(
            tied_weights, tied_values, 2, 1, &maximum_value,
            &selected_indices, &selected_count) ||
        maximum_value != 5 || selected_count != 1 ||
        selected_indices == NULL || selected_indices[0] != 0 ||
        !knapsack_selection_matches(tied_weights, tied_values, 2, 1,
                                    selected_indices, selected_count,
                                    maximum_value, &selected_weight) ||
        selected_weight != 1) {
        goto cleanup;
    }
    free(selected_indices);
    selected_indices = NULL;
    selected_count = 1;
    maximum_value = -1;
    if (dp_knapsack_01_with_selection(
            weights, values, 3, SIZE_MAX, &maximum_value,
            &selected_indices, &selected_count) ||
        maximum_value != 0 || selected_indices != NULL ||
        selected_count != 0) {
        goto cleanup;
    }
    valid = 1;

cleanup:
    free(selected_indices);
    return valid;
}

static int lcs_is_subsequence(const char *candidate, const char *input)
{
    size_t matched = 0;
    size_t i;

    for (i = 0; input[i] != '\0' && candidate[matched] != '\0'; ++i) {
        if (input[i] == candidate[matched]) {
            ++matched;
        }
    }
    return candidate[matched] == '\0';
}

static int run_lcs_case(const char *first,
                        const char *second,
                        size_t expected_length,
                        const char *expected_result)
{
    char *result = NULL;
    size_t length = SIZE_MAX;
    int valid;

    valid = dp_longest_common_subsequence(first, second, &length, &result) &&
            result != NULL && length == expected_length &&
            strlen(result) == length &&
            lcs_is_subsequence(result, first) &&
            lcs_is_subsequence(result, second) &&
            (expected_result == NULL ||
             strcmp(result, expected_result) == 0);
    free(result);
    return valid;
}

static int test_dp_longest_common_subsequence(void)
{
    char first[257];
    char second[257];
    size_t i;

    for (i = 0; i < 256; ++i) {
        first[i] = "ACGT"[i % 4];
        second[i] = first[i];
    }
    first[256] = '\0';
    second[256] = '\0';

    return run_lcs_case("", "", 0, "") &&
           run_lcs_case("", "ACGT", 0, "") &&
           run_lcs_case("GATTACA", "GATTACA", 7, "GATTACA") &&
           run_lcs_case("AAAA", "CCCC", 0, "") &&
           run_lcs_case("ABCBDAB", "BDCABA", 4, NULL) &&
           run_lcs_case("ABC", "BAC", 2, NULL) &&
           run_lcs_case("G", "G", 1, "G") &&
           run_lcs_case("G", "A", 0, "") &&
           run_lcs_case("ACGT", "AGT", 3, "AGT") &&
           run_lcs_case(first, second, 256, first);
}

static int test_dp_longest_common_subsequence_validation(void)
{
    const char embedded_null[] = {'A', '\0', 'B'};
    size_t length = 99;
    char *result = "stale";

    if (dp_longest_common_subsequence(NULL, "A", &length, &result) ||
        length != 0 || result != NULL) {
        return 0;
    }
    length = 99;
    result = "stale";
    if (dp_longest_common_subsequence("A", NULL, &length, &result) ||
        length != 0 || result != NULL) {
        return 0;
    }
    length = 99;
    if (dp_longest_common_subsequence("A", "A", NULL, &result) ||
        result != NULL) {
        return 0;
    }

    length = 99;
    result = "stale";
    if (dp_longest_common_subsequence_n(
            "", 0, embedded_null, sizeof(embedded_null),
            &length, &result) ||
        length != 0 || result != NULL) {
        return 0;
    }
    length = 99;
    result = "stale";
    if (dp_longest_common_subsequence_n(
            embedded_null, sizeof(embedded_null), "", 0,
            &length, &result) ||
        length != 0 || result != NULL) {
        return 0;
    }

    if (!dp_longest_common_subsequence_n(
            "", 0, "ACGT", 4, &length, &result) ||
        length != 0 || result == NULL || result[0] != '\0') {
        free(result);
        return 0;
    }
    free(result);
    result = NULL;
    if (!dp_longest_common_subsequence_n(
            "ACGT", 4, "", 0, &length, &result) ||
        length != 0 || result == NULL || result[0] != '\0') {
        free(result);
        return 0;
    }
    free(result);
    result = NULL;
    if (!dp_longest_common_subsequence_n(
            "", 0, "", 0, &length, &result) ||
        length != 0 || result == NULL || result[0] != '\0') {
        free(result);
        return 0;
    }
    free(result);

    length = 99;
    result = "stale";
    return !dp_longest_common_subsequence_n(
               "A", SIZE_MAX / 2, "B", 3, &length, &result) &&
           length == 0 && result == NULL;
}

static int test_dna_alignment_global(void)
{
    DnaAlignmentResult result = {0};
    int valid;

    valid = dna_alignment_compute("ACGT", "ACGT", DNA_ALIGNMENT_GLOBAL,
                                  2, -1, -2, &result) &&
            result.score == 8 && result.length == 4 &&
            strcmp(result.aligned_first, "ACGT") == 0 &&
            strcmp(result.aligned_second, "ACGT") == 0 &&
            strcmp(result.match_line, "||||") == 0 &&
            result.matches == 4 && result.mismatches == 0 && result.gaps == 0;
    dna_alignment_result_destroy(&result);
    if (!valid) {
        return 0;
    }

    valid = dna_alignment_compute("A", "T", DNA_ALIGNMENT_GLOBAL,
                                  2, -1, -2, &result) &&
            result.score == -1 && result.mismatches == 1 &&
            result.gaps == 0 && strcmp(result.aligned_first, "A") == 0 &&
            strcmp(result.aligned_second, "T") == 0;
    dna_alignment_result_destroy(&result);
    if (!valid) {
        return 0;
    }

    valid = dna_alignment_compute("AAAA", "CCCC", DNA_ALIGNMENT_GLOBAL,
                                  2, -1, -2, &result) &&
            result.score == -4 && result.matches == 0 &&
            result.mismatches == 4 && result.gaps == 0;
    dna_alignment_result_destroy(&result);
    if (!valid) {
        return 0;
    }

    valid = dna_alignment_compute("ACG", "AG", DNA_ALIGNMENT_GLOBAL,
                                  2, -1, -2, &result) &&
            result.score == 2 && result.matches == 2 &&
            result.mismatches == 0 && result.gaps == 1 &&
            strcmp(result.aligned_first, "ACG") == 0 &&
            strcmp(result.aligned_second, "A-G") == 0;
    dna_alignment_result_destroy(&result);
    if (!valid) {
        return 0;
    }

    valid = dna_alignment_compute("ac", "A", DNA_ALIGNMENT_GLOBAL,
                                  2, -1, -2, &result) &&
            result.score == 0 && result.length == 2 &&
            result.matches == 1 && result.gaps == 1;
    dna_alignment_result_destroy(&result);
    return valid;
}

static int test_dna_alignment_local(void)
{
    DnaAlignmentResult result = {0};
    int valid;

    valid = dna_alignment_compute("TTACGTAA", "GGACGCCA", DNA_ALIGNMENT_LOCAL,
                                  2, -1, -2, &result) &&
            result.score == 6 && result.length == 3 &&
            strcmp(result.aligned_first, "ACG") == 0 &&
            strcmp(result.aligned_second, "ACG") == 0 &&
            result.matches == 3 && result.mismatches == 0 &&
            result.gaps == 0 && result.first_start == 2 &&
            result.first_end == 5 && result.second_start == 2 &&
            result.second_end == 5;
    dna_alignment_result_destroy(&result);
    if (!valid) {
        return 0;
    }

    valid = dna_alignment_compute("A", "T", DNA_ALIGNMENT_LOCAL,
                                  2, -1, -2, &result) &&
            result.score == 0 && result.length == 0 &&
            result.aligned_first[0] == '\0' &&
            result.aligned_second[0] == '\0';
    dna_alignment_result_destroy(&result);
    return valid;
}

static int test_dna_alignment_validation_and_overflow(void)
{
    DnaAlignmentResult result = {0};
    int valid;

    if (dna_alignment_compute("", "A", DNA_ALIGNMENT_GLOBAL,
                              2, -1, -2, &result) ||
        dna_alignment_compute("A", "N", DNA_ALIGNMENT_LOCAL,
                              2, -1, -2, &result) ||
        dna_alignment_compute(NULL, "A", DNA_ALIGNMENT_GLOBAL,
                              2, -1, -2, &result)) {
        dna_alignment_result_destroy(&result);
        return 0;
    }
    if (!dna_alignment_compute("AC", "AC", DNA_ALIGNMENT_GLOBAL,
                               2, -1, -2, &result)) {
        return 0;
    }
    valid = !dna_alignment_compute("AA", "AA", DNA_ALIGNMENT_GLOBAL,
                                  INT64_MAX, -1, -2, &result) &&
            result.score == 4 && strcmp(result.aligned_first, "AC") == 0;
    dna_alignment_result_destroy(&result);
    return valid;
}

static int run_matrix_chain_case(const size_t *dimensions,
                                 size_t dimension_count,
                                 uint64_t expected_cost,
                                 const char *expected_order)
{
    uint64_t cost = UINT64_MAX;
    char *order = NULL;
    int valid;

    valid = dp_matrix_chain_multiplication(dimensions, dimension_count,
                                           &cost, &order) &&
            cost == expected_cost && order != NULL &&
            strcmp(order, expected_order) == 0;
    free(order);
    return valid;
}

static int test_dp_matrix_chain_multiplication(void)
{
    const size_t standard[] = {10, 30, 5, 60};
    const size_t single[] = {10, 20};
    const size_t pair[] = {10, 20, 30};
    const size_t tied[] = {1, 1, 1, 1};
    const size_t one_overflowing_order[] = {SIZE_MAX / 2, 1, 3, 1};
    const size_t avoids_overflowing_subchain[] = {
        UINT64_C(2147483648), UINT64_C(2147483648),
        UINT64_C(2147483648), 1, 1
    };
    const size_t all_overflow[] = {SIZE_MAX, SIZE_MAX, SIZE_MAX};
    const size_t zero_dimension[] = {10, 0, 20};
    const size_t one_dimension[] = {10};
    uint64_t cost = UINT64_MAX;
    char *order = "stale";

    if (!run_matrix_chain_case(
            standard, 4, UINT64_C(4500), "((A1 x A2) x A3)") ||
        !run_matrix_chain_case(single, 2, UINT64_C(0), "A1") ||
        !run_matrix_chain_case(pair, 3, UINT64_C(6000), "(A1 x A2)") ||
        !run_matrix_chain_case(
            tied, 4, UINT64_C(2), "(A1 x (A2 x A3))") ||
        !run_matrix_chain_case(
            one_overflowing_order, 4,
            (uint64_t)(SIZE_MAX / 2) + UINT64_C(3),
            "(A1 x (A2 x A3))") ||
        !run_matrix_chain_case(
            avoids_overflowing_subchain, 5,
            (UINT64_C(1) << 63) + UINT64_C(2147483648),
            "(A1 x (A2 x (A3 x A4)))")) {
        return 0;
    }

    if (dp_matrix_chain_multiplication(NULL, 4, &cost, &order) ||
        cost != 0 || order != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    order = "stale";
    if (dp_matrix_chain_multiplication(one_dimension, 1, &cost, &order) ||
        cost != 0 || order != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    order = "stale";
    if (dp_matrix_chain_multiplication(zero_dimension, 3, &cost, &order) ||
        cost != 0 || order != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    order = "stale";
    if (dp_matrix_chain_multiplication(
            standard, SIZE_MAX, &cost, &order) ||
        cost != 0 || order != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    order = "stale";
    return !dp_matrix_chain_multiplication(
               all_overflow, 3, &cost, &order) &&
           cost == 0 && order == NULL;
}

static int obst_test_add(uint64_t left, uint64_t right, uint64_t *sum)
{
    if (left > UINT64_MAX - right) {
        return 0;
    }
    *sum = left + right;
    return 1;
}

static int obst_exhaustive_reference(const uint64_t *frequencies,
                                     size_t first,
                                     size_t last,
                                     uint64_t *minimum_cost)
{
    uint64_t interval_frequency = 0;
    uint64_t best = UINT64_MAX;
    size_t i;
    size_t root;
    int found = 0;

    if (first == last) {
        *minimum_cost = 0;
        return 1;
    }
    for (i = first; i < last; ++i) {
        if (!obst_test_add(interval_frequency, frequencies[i],
                           &interval_frequency)) {
            return 0;
        }
    }
    for (root = first; root < last; ++root) {
        uint64_t left_cost;
        uint64_t right_cost;
        uint64_t child_cost;
        uint64_t candidate;

        if (!obst_exhaustive_reference(frequencies, first, root,
                                       &left_cost) ||
            !obst_exhaustive_reference(frequencies, root + 1, last,
                                       &right_cost) ||
            !obst_test_add(left_cost, right_cost, &child_cost) ||
            !obst_test_add(child_cost, interval_frequency, &candidate)) {
            continue;
        }
        if (!found || candidate < best) {
            best = candidate;
            found = 1;
        }
    }
    if (!found) {
        return 0;
    }
    *minimum_cost = best;
    return 1;
}

static int verify_obst_inorder(const DpOptimalBstNode *node,
                               const int *keys,
                               const uint64_t *frequencies,
                               size_t key_count,
                               size_t *visited,
                               uint64_t depth,
                               uint64_t *search_cost)
{
    uint64_t node_cost;

    if (node == NULL) {
        return 1;
    }
    if (depth == UINT64_MAX && (node->left != NULL || node->right != NULL)) {
        return 0;
    }
    if (!verify_obst_inorder(node->left, keys, frequencies, key_count,
                             visited, depth + 1, search_cost) ||
        *visited >= key_count || node->key != keys[*visited] ||
        node->frequency != frequencies[*visited] ||
        (node->frequency != 0 &&
         depth > UINT64_MAX / node->frequency)) {
        return 0;
    }
    node_cost = node->frequency * depth;
    if (!obst_test_add(*search_cost, node_cost, search_cost)) {
        return 0;
    }
    ++*visited;
    return verify_obst_inorder(node->right, keys, frequencies, key_count,
                               visited, depth + 1, search_cost);
}

static int run_obst_case(const int *keys,
                         const uint64_t *frequencies,
                         size_t key_count,
                         int check_expected,
                         uint64_t expected_cost)
{
    DpOptimalBstNode *root = NULL;
    uint64_t cost = UINT64_MAX;
    uint64_t reference_cost = 0;
    uint64_t reconstructed_cost = 0;
    size_t visited = 0;
    int valid;

    valid = dp_optimal_bst(keys, frequencies, key_count, &cost, &root) &&
            obst_exhaustive_reference(frequencies, 0, key_count,
                                      &reference_cost) &&
            cost == reference_cost &&
            (!check_expected || cost == expected_cost);
    if (valid && key_count != 0) {
        valid = root != NULL &&
                verify_obst_inorder(root, keys, frequencies, key_count,
                                    &visited, 1, &reconstructed_cost) &&
                visited == key_count && reconstructed_cost == cost;
    } else if (valid) {
        valid = root == NULL && cost == 0;
    }
    dp_optimal_bst_destroy(root);
    return valid;
}

static int test_dp_optimal_bst(void)
{
    const int two_keys[] = {10, 20};
    const uint64_t two_frequencies[] = {3, 5};
    const int example_keys[] = {10, 20, 30, 40};
    const uint64_t example_frequencies[] = {34, 8, 50, 3};
    const int zero_keys[] = {2, 4, 6};
    const uint64_t zero_frequencies[] = {0, 5, 0};
    const int all_zero_keys[] = {1, 2, 3, 4};
    const uint64_t all_zero_frequencies[] = {0, 0, 0, 0};
    const int single_key[] = {-7};
    const uint64_t single_frequency[] = {7};
    const int max_frequency_key[] = {99};
    const uint64_t max_frequency[] = {UINT64_MAX};
    const int large_keys[] = {1, 2};
    const uint64_t large_frequencies[] = {
        UINT64_MAX / 4, UINT64_MAX / 8
    };
    const int generated_keys[] = {-20, -10, 0, 10, 20, 30};
    const uint64_t generated_frequencies[] = {3, 0, 7, 2, 5, 1};
    size_t key_count;

    if (!run_obst_case(NULL, NULL, 0, 1, 0) ||
        !run_obst_case(single_key, single_frequency, 1, 1, 7) ||
        !run_obst_case(two_keys, two_frequencies, 2, 1, 11) ||
        !run_obst_case(example_keys, example_frequencies, 4, 1, 148) ||
        !run_obst_case(zero_keys, zero_frequencies, 3, 0, 0) ||
        !run_obst_case(all_zero_keys, all_zero_frequencies, 4, 1, 0) ||
        !run_obst_case(max_frequency_key, max_frequency, 1, 1, UINT64_MAX) ||
        !run_obst_case(large_keys, large_frequencies, 2, 0, 0)) {
        return 0;
    }

    for (key_count = 1; key_count <= 6; ++key_count) {
        if (!run_obst_case(generated_keys, generated_frequencies, key_count,
                           0, 0)) {
            return 0;
        }
    }
    return 1;
}

static int test_dp_optimal_bst_validation(void)
{
    const int unsorted_keys[] = {20, 10};
    const int duplicate_keys[] = {10, 10};
    const int two_keys[] = {1, 2};
    const uint64_t normal_frequencies[] = {1, 2};
    const uint64_t overflowing_prefix[] = {UINT64_MAX, 1};
    const uint64_t overflowing_cost[] = {
        UINT64_MAX / 2, UINT64_MAX / 2
    };
    uint64_t cost = UINT64_MAX;
    DpOptimalBstNode sentinel = {0, 0, NULL, NULL};
    DpOptimalBstNode *root = &sentinel;
    int one_key = 1;
    uint64_t one_frequency = 1;

    if (dp_optimal_bst(unsorted_keys, normal_frequencies, 2, &cost, &root) ||
        cost != 0 || root != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    root = &sentinel;
    if (dp_optimal_bst(duplicate_keys, normal_frequencies, 2, &cost, &root) ||
        cost != 0 || root != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    root = &sentinel;
    if (dp_optimal_bst(NULL, &one_frequency, 1, &cost, &root) ||
        cost != 0 || root != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    root = &sentinel;
    if (dp_optimal_bst(&one_key, NULL, 1, &cost, &root) ||
        cost != 0 || root != NULL) {
        return 0;
    }
    root = &sentinel;
    if (dp_optimal_bst(&one_key, &one_frequency, 1, NULL, &root) ||
        root != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    if (dp_optimal_bst(&one_key, &one_frequency, 1, &cost, NULL) ||
        cost != 0) {
        return 0;
    }
    cost = UINT64_MAX;
    root = &sentinel;
    if (dp_optimal_bst(two_keys, overflowing_prefix, 2, &cost, &root) ||
        cost != 0 || root != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    root = &sentinel;
    if (dp_optimal_bst(two_keys, overflowing_cost, 2, &cost, &root) ||
        cost != 0 || root != NULL) {
        return 0;
    }
    cost = UINT64_MAX;
    root = &sentinel;
    if (dp_optimal_bst(&one_key, &one_frequency, SIZE_MAX, &cost, &root) ||
        cost != 0 || root != NULL) {
        return 0;
    }

    dp_optimal_bst_destroy(NULL);
    return run_obst_case(NULL, NULL, 0, 1, 0);
}

static int bellman_test_add(int64_t left, int64_t right, int64_t *sum)
{
    if ((right > 0 && left > INT64_MAX - right) ||
        (right < 0 && left < INT64_MIN - right)) {
        return 0;
    }
    *sum = left + right;
    return 1;
}

static int bellman_path_matches(const DpBellmanFordEdge *edges,
                                size_t edge_count,
                                const DpBellmanFordResult *result,
                                size_t destination)
{
    size_t *path = NULL;
    size_t path_length = 0;
    size_t i;
    int64_t total_weight = 0;
    int valid;

    if (!dp_bellman_ford_reconstruct_path(
            result, destination, &path, &path_length) ||
        path_length == 0 || path_length > result->vertex_count ||
        path[0] != result->source ||
        path[path_length - 1] != destination) {
        free(path);
        return 0;
    }
    for (i = 0; i + 1 < path_length; ++i) {
        size_t edge_index;
        int found_edge = 0;

        for (edge_index = 0; edge_index < edge_count; ++edge_index) {
            int64_t next_distance;

            if (edges[edge_index].from == path[i] &&
                edges[edge_index].to == path[i + 1] &&
                bellman_test_add(result->distances[path[i]],
                                 edges[edge_index].weight,
                                 &next_distance) &&
                next_distance == result->distances[path[i + 1]] &&
                bellman_test_add(total_weight, edges[edge_index].weight,
                                 &total_weight)) {
                found_edge = 1;
                break;
            }
        }
        if (!found_edge) {
            free(path);
            return 0;
        }
    }
    valid = total_weight == result->distances[destination];
    free(path);
    return valid;
}

static int test_dp_bellman_ford_shortest_paths(void)
{
    const DpBellmanFordEdge known_edges[] = {
        {0, 1, 6}, {0, 2, 7}, {1, 2, 8}, {1, 3, 5}, {1, 4, -4},
        {2, 3, -3}, {2, 4, 9}, {3, 1, -2}, {4, 0, 2}, {4, 3, 7}
    };
    const DpBellmanFordEdge multiple_paths[] = {
        {0, 1, 1}, {0, 2, 1}, {1, 3, 1}, {2, 3, 1}, {3, 4, 0}
    };
    const int64_t known_distances[] = {0, 2, 7, 4, -2};
    DpBellmanFordResult result = {0};
    size_t i;
    int valid;

    valid = dp_bellman_ford(5, known_edges, 10, 0, &result) &&
            !result.negative_cycle_detected &&
            result.vertex_count == 5 && result.source == 0;
    for (i = 0; valid && i < 5; ++i) {
        valid = result.reachable[i] &&
                result.distances[i] == known_distances[i];
    }
    if (valid) {
        valid = bellman_path_matches(known_edges, 10, &result, 4) &&
                result.predecessors[4] == 1;
    }

    if (valid) {
        valid = dp_bellman_ford(6, multiple_paths, 5, 0, &result) &&
                !result.negative_cycle_detected &&
                result.distances[3] == 2 && result.distances[4] == 2 &&
                bellman_path_matches(multiple_paths, 5, &result, 4) &&
                bellman_path_matches(multiple_paths, 5, &result, 0);
    }

    if (valid) {
        valid = !result.reachable[5];
    }

    if (valid) {
        valid = dp_bellman_ford(1, NULL, 0, 0, &result) &&
                !result.negative_cycle_detected &&
                result.reachable[0] && result.distances[0] == 0 &&
                bellman_path_matches(NULL, 0, &result, 0);
    }

    dp_bellman_ford_result_destroy(&result);
    return valid;
}

static int test_dp_bellman_ford_negative_cycles(void)
{
    const DpBellmanFordEdge reachable_cycle[] = {
        {0, 1, 1}, {1, 2, -2}, {2, 1, 1}
    };
    const DpBellmanFordEdge unreachable_cycle[] = {
        {0, 1, 2}, {2, 3, -2}, {3, 2, 1}
    };
    DpBellmanFordResult result = {0};
    size_t *path = NULL;
    size_t path_length = 0;
    int valid;

    valid = dp_bellman_ford(3, reachable_cycle, 3, 0, &result) &&
            result.negative_cycle_detected &&
            !dp_bellman_ford_reconstruct_path(
                &result, 0, &path, &path_length) &&
            path == NULL && path_length == 0;
    if (valid) {
        valid = dp_bellman_ford(4, unreachable_cycle, 3, 0, &result) &&
                !result.negative_cycle_detected &&
                result.reachable[0] && result.reachable[1] &&
                !result.reachable[2] && !result.reachable[3] &&
                result.distances[1] == 2 &&
                bellman_path_matches(unreachable_cycle, 3, &result, 1);
    }
    dp_bellman_ford_result_destroy(&result);
    return valid;
}

static int test_dp_bellman_ford_validation_and_reuse(void)
{
    const DpBellmanFordEdge basic_edges[] = {{0, 1, 4}, {1, 2, -1}};
    const DpBellmanFordEdge invalid_edge[] = {{0, 3, 1}};
    const DpBellmanFordEdge overflowing_distance[] = {
        {0, 1, INT64_MAX}, {1, 2, 1}
    };
    const DpBellmanFordEdge overflowing_negative[] = {
        {0, 1, INT64_MIN}, {1, 2, -1}
    };
    DpBellmanFordResult result = {0};
    int64_t *saved_distances;
    size_t *saved_predecessors;
    unsigned char *saved_reachable;
    size_t *path = NULL;
    size_t path_length = 99;
    size_t i;
    int valid;

    if (!dp_bellman_ford(3, basic_edges, 2, 0, &result)) {
        return 0;
    }
    saved_distances = result.distances;
    saved_predecessors = result.predecessors;
    saved_reachable = result.reachable;

    if (dp_bellman_ford(3, basic_edges, 2, 3, &result) ||
        dp_bellman_ford(3, invalid_edge, 1, 0, &result) ||
        dp_bellman_ford(0, NULL, 0, 0, &result) ||
        dp_bellman_ford(3, NULL, 1, 0, &result) ||
        dp_bellman_ford(3, basic_edges, 2, 0, NULL) ||
        result.distances != saved_distances ||
        result.predecessors != saved_predecessors ||
        result.reachable != saved_reachable ||
        result.distances[2] != 3) {
        dp_bellman_ford_result_destroy(&result);
        return 0;
    }

    if (dp_bellman_ford_reconstruct_path(
            &result, 3, &path, &path_length) ||
        path != NULL || path_length != 0 ||
        dp_bellman_ford_reconstruct_path(
            &result, 2, NULL, &path_length) ||
        path_length != 0) {
        dp_bellman_ford_result_destroy(&result);
        return 0;
    }

    result.predecessors[1] = 2;
    result.predecessors[2] = 1;
    if (dp_bellman_ford_reconstruct_path(
            &result, 2, &path, &path_length) ||
        path != NULL || path_length != 0) {
        dp_bellman_ford_result_destroy(&result);
        return 0;
    }
    result.predecessors[1] = 0;
    result.predecessors[2] = 1;

    if (dp_bellman_ford(3, overflowing_distance, 2, 0, &result) ||
        dp_bellman_ford(3, overflowing_negative, 2, 0, &result) ||
        dp_bellman_ford(SIZE_MAX, NULL, 0, 0, &result) ||
        dp_bellman_ford(2, basic_edges,
                        SIZE_MAX / sizeof(*basic_edges) + 1, 0, &result) ||
        result.distances != saved_distances ||
        result.predecessors != saved_predecessors ||
        result.reachable != saved_reachable) {
        dp_bellman_ford_result_destroy(&result);
        return 0;
    }

    valid = 1;
    for (i = 0; i < 16; ++i) {
        valid = dp_bellman_ford(3, basic_edges, 2, 0, &result) &&
                result.distances[2] == 3 &&
                !result.negative_cycle_detected;
        if (!valid) {
            break;
        }
    }
    dp_bellman_ford_result_destroy(&result);
    return valid && result.vertex_count == 0 && result.distances == NULL &&
           result.predecessors == NULL && result.reachable == NULL;
}

static int floyd_path_matches(const DpWeightedDirectedMatrix *graph,
                              const DpFloydWarshallResult *result,
                              size_t source,
                              size_t destination)
{
    size_t *path = NULL;
    size_t path_length = 0;
    size_t i;
    int64_t total = 0;
    int valid;

    if (!dp_floyd_warshall_reconstruct_path(
            result, source, destination, &path, &path_length) ||
        path_length == 0 || path_length > result->vertex_count ||
        path[0] != source || path[path_length - 1] != destination) {
        free(path);
        return 0;
    }
    for (i = 0; i + 1 < path_length; ++i) {
        size_t cell = path[i] * graph->vertex_count + path[i + 1];

        if (!graph->edge_present[cell] ||
            !bellman_test_add(total, graph->weights[cell], &total)) {
            free(path);
            return 0;
        }
    }
    valid = total == result->distances[
                         source * result->vertex_count + destination];
    free(path);
    return valid;
}

static int test_dp_floyd_warshall_paths(void)
{
    const int64_t single_weights[] = {0};
    const unsigned char single_present[] = {0};
    const int64_t two_weights[] = {0, 3, 0, 0};
    const unsigned char two_present[] = {0, 1, 0, 0};
    const int64_t multiple_weights[] = {
        0, 5, 10, 0,
        0, 0, 0, 7,
        0, 0, 0, -2,
        0, 0, 0, 0
    };
    const unsigned char multiple_present[] = {
        0, 1, 1, 0,
        0, 0, 1, 1,
        0, 0, 0, 1,
        0, 0, 0, 0
    };
    DpWeightedDirectedMatrix single = {1, single_weights, single_present};
    DpWeightedDirectedMatrix two = {2, two_weights, two_present};
    DpWeightedDirectedMatrix multiple = {
        4, multiple_weights, multiple_present
    };
    DpFloydWarshallResult result = {0};
    size_t *path = NULL;
    size_t path_length = 0;
    int valid;

    valid = dp_floyd_warshall(&single, &result) &&
            result.distances[0] == 0 && result.reachable[0] &&
            !result.negative_cycle_detected &&
            floyd_path_matches(&single, &result, 0, 0);
    if (valid) {
        valid = dp_floyd_warshall(&two, &result) &&
                result.distances[1] == 3 &&
                result.reachable[1] && !result.reachable[2] &&
                floyd_path_matches(&two, &result, 0, 1) &&
                !dp_floyd_warshall_reconstruct_path(
                    &result, 1, 0, &path, &path_length) &&
                path == NULL && path_length == 0;
    }
    if (valid) {
        valid = dp_floyd_warshall(&multiple, &result) &&
                result.distances[2] == 5 && result.distances[3] == 3 &&
                result.reachable[3] &&
                floyd_path_matches(&multiple, &result, 0, 3) &&
                floyd_path_matches(&multiple, &result, 2, 3) &&
                floyd_path_matches(&multiple, &result, 3, 3);
    }
    dp_floyd_warshall_result_destroy(&result);
    return valid;
}

static int test_dp_floyd_warshall_negative_cycles(void)
{
    const int64_t reachable_weights[] = {
        0, 1, 0,
        0, 0, -2,
        0, 1, 0
    };
    const unsigned char reachable_present[] = {
        0, 1, 0,
        0, 0, 1,
        0, 1, 0
    };
    const int64_t unreachable_weights[] = {
        0, 2, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, -2,
        0, 0, 1, 0
    };
    const unsigned char unreachable_present[] = {
        0, 1, 0, 0,
        0, 0, 0, 0,
        0, 0, 0, 1,
        0, 0, 1, 0
    };
    DpWeightedDirectedMatrix reachable_cycle = {
        3, reachable_weights, reachable_present
    };
    DpWeightedDirectedMatrix unreachable_cycle = {
        4, unreachable_weights, unreachable_present
    };
    DpFloydWarshallResult result = {0};
    size_t *path = NULL;
    size_t path_length = 0;
    size_t affected_cell;
    int valid;

    valid = dp_floyd_warshall(&reachable_cycle, &result) &&
            result.negative_cycle_detected;
    affected_cell = 0 * result.vertex_count + 2;
    if (valid) {
        valid = result.negative_cycle_affected[affected_cell] &&
                !dp_floyd_warshall_reconstruct_path(
                    &result, 0, 2, &path, &path_length) &&
                path == NULL && path_length == 0 &&
                floyd_path_matches(&reachable_cycle, &result, 0, 0);
    }
    if (valid) {
        valid = dp_floyd_warshall(&unreachable_cycle, &result) &&
                result.negative_cycle_detected &&
                result.reachable[1] && !result.reachable[2] &&
                result.negative_cycle_affected[2 * 4 + 3] &&
                !result.negative_cycle_affected[0 * 4 + 1] &&
                floyd_path_matches(&unreachable_cycle, &result, 0, 1) &&
                !dp_floyd_warshall_reconstruct_path(
                    &result, 0, 2, &path, &path_length);
    }
    dp_floyd_warshall_result_destroy(&result);
    return valid;
}

static int test_dp_floyd_warshall_validation_and_reuse(void)
{
    const int64_t valid_weights[] = {
        0, 2, 0,
        0, 0, 1,
        0, 0, 0
    };
    const unsigned char valid_present[] = {
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    };
    const int64_t overflow_weights[] = {
        0, INT64_MAX, 0,
        0, 0, 1,
        0, 0, 0
    };
    const unsigned char overflow_present[] = {
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    };
    DpWeightedDirectedMatrix valid_graph = {
        3, valid_weights, valid_present
    };
    DpWeightedDirectedMatrix overflow_graph = {
        3, overflow_weights, overflow_present
    };
    DpWeightedDirectedMatrix null_matrix = {3, NULL, valid_present};
    DpWeightedDirectedMatrix oversized = {SIZE_MAX, valid_weights,
                                          valid_present};
    DpFloydWarshallResult result = {0};
    int64_t *saved_distances;
    size_t *saved_next_hops;
    unsigned char *saved_reachable;
    unsigned char *saved_affected;
    size_t *path = NULL;
    size_t path_length = 42;
    size_t i;
    int valid;

    if (!dp_floyd_warshall(&valid_graph, &result)) {
        return 0;
    }
    saved_distances = result.distances;
    saved_next_hops = result.next_hops;
    saved_reachable = result.reachable;
    saved_affected = result.negative_cycle_affected;

    if (dp_floyd_warshall(NULL, &result) ||
        dp_floyd_warshall(&null_matrix, &result) ||
        dp_floyd_warshall(&oversized, &result) ||
        dp_floyd_warshall(&valid_graph, NULL) ||
        dp_floyd_warshall(&overflow_graph, &result) ||
        result.distances != saved_distances ||
        result.next_hops != saved_next_hops ||
        result.reachable != saved_reachable ||
        result.negative_cycle_affected != saved_affected) {
        dp_floyd_warshall_result_destroy(&result);
        return 0;
    }

    if (dp_floyd_warshall_reconstruct_path(
            &result, 3, 0, &path, &path_length) ||
        path != NULL || path_length != 0 ||
        dp_floyd_warshall_reconstruct_path(
            &result, 0, 3, &path, &path_length) ||
        path != NULL || path_length != 0 ||
        dp_floyd_warshall_reconstruct_path(
            &result, 0, 2, NULL, &path_length) ||
        path_length != 0) {
        dp_floyd_warshall_result_destroy(&result);
        return 0;
    }

    result.next_hops[0 * 3 + 2] = 1;
    result.next_hops[1 * 3 + 2] = 0;
    if (dp_floyd_warshall_reconstruct_path(
            &result, 0, 2, &path, &path_length) ||
        path != NULL || path_length != 0) {
        dp_floyd_warshall_result_destroy(&result);
        return 0;
    }

    valid = 1;
    for (i = 0; i < 8; ++i) {
        valid = dp_floyd_warshall(&valid_graph, &result) &&
                result.distances[2] == 3 &&
                !result.negative_cycle_detected;
        if (!valid) {
            break;
        }
    }
    dp_floyd_warshall_result_destroy(&result);
    return valid && result.vertex_count == 0 && result.distances == NULL &&
           result.next_hops == NULL && result.reachable == NULL &&
           result.negative_cycle_affected == NULL;
}

static int compare_bellman_ford_with_floyd(
    size_t vertex_count,
    const DpBellmanFordEdge *edges,
    size_t edge_count)
{
    int64_t weights[36] = {0};
    unsigned char present[36] = {0};
    DpWeightedDirectedMatrix graph;
    DpFloydWarshallResult all_pairs = {0};
    DpBellmanFordResult single_source = {0};
    size_t i;
    size_t source;
    int valid = 0;

    if (vertex_count == 0 || vertex_count > 6) {
        return 0;
    }
    for (i = 0; i < edge_count; ++i) {
        size_t cell = edges[i].from * vertex_count + edges[i].to;

        if (!present[cell] || edges[i].weight < weights[cell]) {
            weights[cell] = edges[i].weight;
            present[cell] = 1;
        }
    }
    graph.vertex_count = vertex_count;
    graph.weights = weights;
    graph.edge_present = present;
    if (!dp_floyd_warshall(&graph, &all_pairs) ||
        all_pairs.negative_cycle_detected) {
        goto cleanup;
    }

    valid = 1;
    for (source = 0; valid && source < vertex_count; ++source) {
        size_t destination;

        if (!dp_bellman_ford(vertex_count, edges, edge_count, source,
                             &single_source) ||
            single_source.negative_cycle_detected) {
            valid = 0;
            break;
        }
        for (destination = 0; destination < vertex_count; ++destination) {
            size_t cell = source * vertex_count + destination;

            if (single_source.reachable[destination] !=
                    all_pairs.reachable[cell] ||
                (all_pairs.reachable[cell] &&
                 single_source.distances[destination] !=
                     all_pairs.distances[cell])) {
                valid = 0;
                break;
            }
        }
    }

cleanup:
    dp_bellman_ford_result_destroy(&single_source);
    dp_floyd_warshall_result_destroy(&all_pairs);
    return valid;
}

static int test_dp_floyd_bellman_ford_agreement(void)
{
    const DpBellmanFordEdge positive_edges[] = {
        {0, 1, 4}, {0, 2, 10}, {1, 2, 0}, {2, 3, 3}, {1, 3, 8}
    };
    const DpBellmanFordEdge negative_edges[] = {
        {0, 1, 6}, {0, 2, 7}, {1, 2, 8}, {1, 3, 5}, {1, 4, -4},
        {2, 3, -3}, {2, 4, 9}, {3, 1, -2}, {4, 0, 2}, {4, 3, 7}
    };
    const DpBellmanFordEdge disconnected_edges[] = {
        {0, 1, 0}, {1, 2, -2}, {3, 4, 5}
    };

    return compare_bellman_ford_with_floyd(
               4, positive_edges, 5) &&
           compare_bellman_ford_with_floyd(
               5, negative_edges, 10) &&
           compare_bellman_ford_with_floyd(
               6, disconnected_edges, 3);
}

static int tsp_brute_force_visit(const DpWeightedDirectedMatrix *graph,
                                 size_t start,
                                 size_t last,
                                 uint64_t visited,
                                 size_t visited_count,
                                 int64_t cost,
                                 int64_t *best_cost,
                                 int *found_tour)
{
    size_t vertex;

    if (visited_count == graph->vertex_count) {
        size_t closing_edge = last * graph->vertex_count + start;
        int64_t completed_cost;

        if (graph->edge_present[closing_edge] &&
            bellman_test_add(cost, graph->weights[closing_edge],
                             &completed_cost) &&
            (!*found_tour || completed_cost < *best_cost)) {
            *best_cost = completed_cost;
            *found_tour = 1;
        }
        return 1;
    }

    for (vertex = 0; vertex < graph->vertex_count; ++vertex) {
        uint64_t bit = UINT64_C(1) << vertex;
        size_t edge_cell = last * graph->vertex_count + vertex;
        int64_t next_cost;

        if ((visited & bit) != 0 || !graph->edge_present[edge_cell] ||
            !bellman_test_add(cost, graph->weights[edge_cell], &next_cost)) {
            continue;
        }
        if (!tsp_brute_force_visit(graph, start, vertex, visited | bit,
                                   visited_count + 1, next_cost,
                                   best_cost, found_tour)) {
            return 0;
        }
    }
    return 1;
}

static int tsp_brute_force_reference(const DpWeightedDirectedMatrix *graph,
                                     size_t start,
                                     int64_t *minimum_cost,
                                     int *found_tour)
{
    uint64_t start_bit = UINT64_C(1) << start;

    *found_tour = 0;
    *minimum_cost = 0;
    return tsp_brute_force_visit(graph, start, start, start_bit, 1, 0,
                                 minimum_cost, found_tour);
}

static int tsp_result_is_valid(const DpWeightedDirectedMatrix *graph,
                               size_t start,
                               const DpTspResult *result)
{
    unsigned char seen[DP_TSP_MAX_VERTICES] = {0};
    int64_t total = 0;
    size_t i;

    if (result->vertex_count != graph->vertex_count ||
        result->tour == NULL ||
        result->tour_length != graph->vertex_count + 1 ||
        result->tour[0] != start ||
        result->tour[result->tour_length - 1] != start) {
        return 0;
    }
    for (i = 0; i < graph->vertex_count; ++i) {
        size_t vertex = result->tour[i];

        if (vertex >= graph->vertex_count || seen[vertex]) {
            return 0;
        }
        seen[vertex] = 1;
    }
    for (i = 0; i + 1 < result->tour_length; ++i) {
        size_t cell = result->tour[i] * graph->vertex_count +
                      result->tour[i + 1];

        if (!graph->edge_present[cell] ||
            !bellman_test_add(total, graph->weights[cell], &total)) {
            return 0;
        }
    }
    return total == result->minimum_cost;
}

static int run_tsp_reference_case(const DpWeightedDirectedMatrix *graph,
                                  size_t start,
                                  int64_t known_cost)
{
    DpTspResult result = {0};
    int64_t reference_cost;
    int found_tour;
    DpTspStatus status;
    int valid;

    if (!tsp_brute_force_reference(graph, start, &reference_cost,
                                   &found_tour)) {
        return 0;
    }
    status = dp_tsp_bitmask(graph, start, &result);
    if (!found_tour) {
        valid = status == DP_TSP_NO_TOUR &&
                result.tour == NULL && result.tour_length == 0 &&
                result.vertex_count == 0;
    } else {
        valid = status == DP_TSP_TOUR_FOUND &&
                result.minimum_cost == reference_cost &&
                result.minimum_cost == known_cost &&
                tsp_result_is_valid(graph, start, &result);
    }
    dp_tsp_result_destroy(&result);
    return valid;
}

static int test_dp_tsp_bitmask_examples(void)
{
    const int64_t single_weights[] = {9};
    const unsigned char single_present[] = {1};
    const int64_t two_weights[] = {0, 3, 5, 0};
    const unsigned char two_present[] = {0, 1, 1, 0};
    const int64_t standard_weights[] = {
        0, 10, 15, 20,
        5, 0, 9, 10,
        6, 13, 0, 12,
        8, 8, 9, 0
    };
    const unsigned char standard_present[] = {
        0, 1, 1, 1,
        1, 0, 1, 1,
        1, 1, 0, 1,
        1, 1, 1, 0
    };
    const int64_t negative_weights[] = {
        0, -2, 10,
        3, 0, 4,
        1, -5, 0
    };
    const unsigned char negative_present[] = {
        0, 1, 1,
        1, 0, 1,
        1, 1, 0
    };
    DpWeightedDirectedMatrix single = {
        1, single_weights, single_present
    };
    DpWeightedDirectedMatrix two = {2, two_weights, two_present};
    DpWeightedDirectedMatrix standard = {
        4, standard_weights, standard_present
    };
    DpWeightedDirectedMatrix negative = {
        3, negative_weights, negative_present
    };

    return run_tsp_reference_case(&single, 0, 9) &&
           run_tsp_reference_case(&two, 0, 8) &&
           run_tsp_reference_case(&standard, 0, 35) &&
           run_tsp_reference_case(&negative, 0, 3);
}

static int test_dp_tsp_no_tour_and_limits(void)
{
    const int64_t partial_weights[] = {
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    };
    const unsigned char partial_present[] = {
        0, 1, 0,
        0, 0, 1,
        0, 0, 0
    };
    int64_t maximum_weights[
        DP_TSP_MAX_VERTICES * DP_TSP_MAX_VERTICES] = {0};
    unsigned char maximum_present[
        DP_TSP_MAX_VERTICES * DP_TSP_MAX_VERTICES] = {0};
    DpWeightedDirectedMatrix partial = {
        3, partial_weights, partial_present
    };
    DpWeightedDirectedMatrix maximum = {
        DP_TSP_MAX_VERTICES, maximum_weights, maximum_present
    };
    DpTspResult result = {0};
    DpTspStatus status;

    if (!run_tsp_reference_case(&partial, 0, 0)) {
        return 0;
    }
    status = dp_tsp_bitmask(&maximum, 0, &result);
    dp_tsp_result_destroy(&result);
    return status == DP_TSP_NO_TOUR && result.tour == NULL &&
           result.tour_length == 0 && result.vertex_count == 0;
}

static int test_dp_tsp_validation_and_reuse(void)
{
    const int64_t valid_weights[] = {0, 2, 5, 0};
    const unsigned char valid_present[] = {0, 1, 1, 0};
    const int64_t overflow_weights[] = {
        0, INT64_MAX,
        1, 0
    };
    const unsigned char overflow_present[] = {0, 1, 1, 0};
    DpWeightedDirectedMatrix valid_graph = {
        2, valid_weights, valid_present
    };
    DpWeightedDirectedMatrix overflow_graph = {
        2, overflow_weights, overflow_present
    };
    DpWeightedDirectedMatrix null_matrix = {2, NULL, valid_present};
    DpWeightedDirectedMatrix too_large = {
        DP_TSP_MAX_VERTICES + 1, valid_weights, valid_present
    };
    DpWeightedDirectedMatrix oversized = {
        SIZE_MAX, valid_weights, valid_present
    };
    DpTspResult result = {0};
    size_t *saved_tour;
    size_t i;
    int valid;

    if (dp_tsp_bitmask(NULL, 0, &result) != DP_TSP_ERROR ||
        dp_tsp_bitmask(&valid_graph, 2, &result) != DP_TSP_ERROR ||
        dp_tsp_bitmask(&null_matrix, 0, &result) != DP_TSP_ERROR ||
        dp_tsp_bitmask(&too_large, 0, &result) != DP_TSP_ERROR ||
        dp_tsp_bitmask(&oversized, 0, &result) != DP_TSP_ERROR ||
        dp_tsp_bitmask(&valid_graph, 0, NULL) != DP_TSP_ERROR ||
        dp_tsp_bitmask(&valid_graph, 0, &result) != DP_TSP_TOUR_FOUND) {
        dp_tsp_result_destroy(&result);
        return 0;
    }
    saved_tour = result.tour;
    if (dp_tsp_bitmask(&overflow_graph, 0, &result) != DP_TSP_ERROR ||
        result.tour != saved_tour || result.minimum_cost != 7 ||
        dp_tsp_bitmask(&valid_graph, 1, &result) != DP_TSP_TOUR_FOUND ||
        result.tour[0] != 1 || result.tour[result.tour_length - 1] != 1) {
        dp_tsp_result_destroy(&result);
        return 0;
    }

    valid = 1;
    for (i = 0; i < 8; ++i) {
        valid = dp_tsp_bitmask(&valid_graph, 0, &result) ==
                    DP_TSP_TOUR_FOUND &&
                tsp_result_is_valid(&valid_graph, 0, &result);
        if (!valid) {
            break;
        }
    }
    dp_tsp_result_destroy(&result);
    return valid && result.tour == NULL && result.tour_length == 0 &&
           result.vertex_count == 0;
}

static int sum_of_subsets_result_is_valid(const int64_t *values,
                                         size_t item_count,
                                         int64_t target,
                                         const SumOfSubsetsResult *result)
{
    size_t subset_index;

    for (subset_index = 0; subset_index < result->subset_count;
         ++subset_index) {
        const SumOfSubsetsSubset *subset = &result->subsets[subset_index];
        int64_t sum = 0;
        size_t i;
        size_t earlier;

        if (subset->index_count > item_count ||
            ((subset->index_count == 0) != (subset->indices == NULL))) {
            return 0;
        }
        for (i = 0; i < subset->index_count; ++i) {
            if (subset->indices[i] >= item_count ||
                (i != 0 && subset->indices[i - 1] >= subset->indices[i]) ||
                values[subset->indices[i]] <= 0 ||
                values[subset->indices[i]] > target - sum) {
                return 0;
            }
            sum += values[subset->indices[i]];
        }
        if (sum != target) {
            return 0;
        }
        for (earlier = 0; earlier < subset_index; ++earlier) {
            const SumOfSubsetsSubset *other = &result->subsets[earlier];

            if (other->index_count == subset->index_count &&
                (subset->index_count == 0 ||
                 memcmp(other->indices, subset->indices,
                        subset->index_count * sizeof(*subset->indices)) == 0)) {
                return 0;
            }
        }
    }
    return 1;
}

static int sum_subset_contains(const SumOfSubsetsResult *result,
                               const size_t *indices,
                               size_t index_count)
{
    size_t i;

    for (i = 0; i < result->subset_count; ++i) {
        if (result->subsets[i].index_count == index_count &&
            (index_count == 0 ||
             memcmp(result->subsets[i].indices, indices,
                    index_count * sizeof(*indices)) == 0)) {
            return 1;
        }
    }
    return 0;
}

static int test_backtracking_sum_of_subsets(void)
{
    const int64_t standard[] = {5, 10, 12, 13, 15, 18};
    const int64_t duplicates[] = {2, 2, 3};
    const int64_t single[] = {7};
    const int64_t invalid[] = {2, -1, 4};
    const int64_t zero_input[] = {2, 0, 3};
    const size_t first_duplicate_solution[] = {0, 2};
    const size_t second_duplicate_solution[] = {1, 2};
    const size_t pair_solution[] = {0, 1};
    SumOfSubsetsResult result = {NULL, 0};
    SumOfSubsetsSubset *saved_subsets;
    int valid;

    if (!backtracking_sum_of_subsets(standard, 6, 30, &result) ||
        result.subset_count != 3 ||
        !sum_of_subsets_result_is_valid(standard, 6, 30, &result) ||
        !backtracking_sum_of_subsets(duplicates, 3, 5, &result) ||
        result.subset_count != 2 ||
        !sum_of_subsets_result_is_valid(duplicates, 3, 5, &result) ||
        !sum_subset_contains(&result, first_duplicate_solution, 2) ||
        !sum_subset_contains(&result, second_duplicate_solution, 2) ||
        !backtracking_sum_of_subsets(duplicates, 3, 4, &result) ||
        result.subset_count != 1 ||
        !sum_subset_contains(&result, pair_solution, 2) ||
        !backtracking_sum_of_subsets(single, 1, 7, &result) ||
        result.subset_count != 1 ||
        !backtracking_sum_of_subsets(single, 1, 8, &result) ||
        result.subset_count != 0 ||
        !backtracking_sum_of_subsets(NULL, 0, 1, &result) ||
        result.subset_count != 0 ||
        !backtracking_sum_of_subsets(NULL, 0, 0, &result) ||
        result.subset_count != 1 ||
        result.subsets[0].indices != NULL ||
        result.subsets[0].index_count != 0 ||
        !backtracking_sum_of_subsets(standard, 6, 0, &result) ||
        result.subset_count != 1 ||
        result.subsets[0].index_count != 0) {
        backtracking_sum_of_subsets_result_destroy(&result);
        return 0;
    }

    saved_subsets = result.subsets;
    if (backtracking_sum_of_subsets(invalid, 3, 4, &result) ||
        backtracking_sum_of_subsets(zero_input, 3, 4, &result) ||
        backtracking_sum_of_subsets(standard, 6, -1, &result) ||
        backtracking_sum_of_subsets(NULL, 1, 2, &result) ||
        backtracking_sum_of_subsets(standard, SIZE_MAX, 2, &result) ||
        backtracking_sum_of_subsets(standard, 6, 4, NULL) ||
        result.subsets != saved_subsets || result.subset_count != 1) {
        backtracking_sum_of_subsets_result_destroy(&result);
        return 0;
    }
    valid = sum_of_subsets_result_is_valid(standard, 6, 0, &result);
    backtracking_sum_of_subsets_result_destroy(&result);
    return valid && result.subsets == NULL && result.subset_count == 0;
}

typedef int (*SearchKnapsackFunction)(const size_t *,
                                      const int64_t *,
                                      size_t,
                                      size_t,
                                      KnapsackSearchResult *);

static int search_knapsack_result_is_valid(const size_t *weights,
                                           const int64_t *values,
                                           size_t item_count,
                                           size_t capacity,
                                           const KnapsackSearchResult *result)
{
    size_t total_weight = 0;
    int64_t total_value = 0;
    size_t i;

    if (result->selected_count > item_count ||
        ((result->selected_count == 0) !=
         (result->selected_indices == NULL))) {
        return 0;
    }
    for (i = 0; i < result->selected_count; ++i) {
        size_t item = result->selected_indices[i];

        if (item >= item_count ||
            (i != 0 && result->selected_indices[i - 1] >= item) ||
            weights[item] > capacity - total_weight ||
            values[item] > INT64_MAX - total_value) {
            return 0;
        }
        total_weight += weights[item];
        total_value += values[item];
    }
    return total_weight <= capacity &&
           total_value == result->maximum_value;
}

static int run_search_knapsack_case(SearchKnapsackFunction solve,
                                    const size_t *weights,
                                    const int64_t *values,
                                    size_t item_count,
                                    size_t capacity,
                                    int64_t expected_value)
{
    KnapsackSearchResult result = {0, NULL, 0};
    int valid = solve(weights, values, item_count, capacity, &result) &&
                result.maximum_value == expected_value &&
                search_knapsack_result_is_valid(weights, values, item_count,
                                                capacity, &result);

    knapsack_search_result_destroy(&result);
    return valid;
}

static int run_search_knapsack_common_cases(SearchKnapsackFunction solve)
{
    const size_t standard_weights[] = {10, 20, 30};
    const int64_t standard_values[] = {60, 100, 120};
    const size_t zero_weights[] = {0, 0, 2, 4};
    const int64_t zero_values[] = {7, 8, 2, 9};
    const size_t no_fit_weights[] = {3, 5};
    const int64_t no_fit_values[] = {4, 7};
    const size_t tie_weights[] = {2, 2};
    const int64_t tie_values[] = {3, 3};
    const size_t all_fit_weights[] = {1, 2, 0};
    const int64_t all_fit_values[] = {5, 6, 4};

    return run_search_knapsack_case(solve, standard_weights, standard_values,
                                    3, 50, 220) &&
           run_search_knapsack_case(solve, zero_weights, zero_values,
                                    4, 4, 24) &&
           run_search_knapsack_case(solve, no_fit_weights, no_fit_values,
                                    2, 0, 0) &&
           run_search_knapsack_case(solve, tie_weights, tie_values, 2, 2, 3) &&
           run_search_knapsack_case(solve, all_fit_weights, all_fit_values,
                                    3, 3, 15) &&
           run_search_knapsack_case(solve, NULL, NULL, 0, SIZE_MAX, 0);
}

static int run_search_knapsack_failures(SearchKnapsackFunction solve)
{
    const size_t weights[] = {0, 0};
    const int64_t valid_values[] = {3, 4};
    const int64_t invalid_values[] = {3, -1};
    const int64_t overflow_values[] = {INT64_MAX, 1};
    KnapsackSearchResult result = {0, NULL, 0};
    size_t *saved_indices;
    int valid;

    if (!solve(weights, valid_values, 2, 0, &result)) {
        return 0;
    }
    saved_indices = result.selected_indices;
    if (solve(NULL, valid_values, 2, 0, &result) ||
        solve(weights, invalid_values, 2, 0, &result) ||
        solve(weights, overflow_values, 2, 0, &result) ||
        solve(weights, valid_values, SIZE_MAX, 0, &result) ||
        solve(weights, valid_values, 2, 0, NULL) ||
        result.selected_indices != saved_indices ||
        result.maximum_value != 7 || result.selected_count != 2) {
        knapsack_search_result_destroy(&result);
        return 0;
    }
    valid = search_knapsack_result_is_valid(
        weights, valid_values, 2, 0, &result);
    knapsack_search_result_destroy(&result);
    return valid;
}

static int test_backtracking_knapsack(void)
{
    return run_search_knapsack_common_cases(backtracking_knapsack_01) &&
           run_search_knapsack_failures(backtracking_knapsack_01);
}

static int test_fifo_branch_bound_knapsack(void)
{
    return run_search_knapsack_common_cases(branch_bound_fifo_knapsack_01) &&
           run_search_knapsack_failures(branch_bound_fifo_knapsack_01);
}

static int test_lifo_branch_bound_knapsack(void)
{
    return run_search_knapsack_common_cases(branch_bound_lifo_knapsack_01) &&
           run_search_knapsack_failures(branch_bound_lifo_knapsack_01);
}

static int test_least_cost_branch_bound_knapsack(void)
{
    return run_search_knapsack_common_cases(
               branch_bound_least_cost_knapsack_01) &&
           run_search_knapsack_failures(
               branch_bound_least_cost_knapsack_01);
}

static int64_t brute_force_knapsack_value(const size_t *weights,
                                         const int64_t *values,
                                         size_t item_count,
                                         size_t capacity)
{
    uint64_t subset_count = UINT64_C(1) << item_count;
    int64_t best_value = 0;
    uint64_t mask;

    for (mask = 0; mask < subset_count; ++mask) {
        size_t total_weight = 0;
        int64_t total_value = 0;
        size_t item;
        int feasible = 1;

        for (item = 0; item < item_count; ++item) {
            if ((mask & (UINT64_C(1) << item)) == 0) {
                continue;
            }
            if (weights[item] > capacity - total_weight ||
                values[item] > INT64_MAX - total_value) {
                feasible = 0;
                break;
            }
            total_weight += weights[item];
            total_value += values[item];
        }
        if (feasible && total_value > best_value) {
            best_value = total_value;
        }
    }
    return best_value;
}

static uint32_t next_knapsack_test_random(uint32_t *state)
{
    *state = *state * UINT32_C(1664525) + UINT32_C(1013904223);
    return *state;
}

static int test_knapsack_algorithm_agreement(void)
{
    uint32_t random_state = UINT32_C(0x5eed1234);
    size_t test_case;

    for (test_case = 0; test_case < 64; ++test_case) {
        size_t weights[8];
        int64_t values[8];
        size_t item_count =
            (size_t)(next_knapsack_test_random(&random_state) % 9);
        size_t capacity =
            (size_t)(next_knapsack_test_random(&random_state) % 31);
        int64_t dp_value;
        int64_t expected;
        KnapsackSearchResult results[4] = {
            {0, NULL, 0}, {0, NULL, 0}, {0, NULL, 0}, {0, NULL, 0}
        };
        size_t i;
        int valid = 1;
        SearchKnapsackFunction solvers[] = {
            backtracking_knapsack_01,
            branch_bound_fifo_knapsack_01,
            branch_bound_lifo_knapsack_01,
            branch_bound_least_cost_knapsack_01
        };

        for (i = 0; i < item_count; ++i) {
            weights[i] =
                (size_t)(next_knapsack_test_random(&random_state) % 12);
            values[i] =
                (int64_t)(next_knapsack_test_random(&random_state) % 25);
        }
        expected = brute_force_knapsack_value(
            weights, values, item_count, capacity);
        if (!dp_knapsack_01(weights, values, item_count, capacity, &dp_value) ||
            dp_value != expected) {
            valid = 0;
        }
        for (i = 0; valid && i < 4; ++i) {
            valid = solvers[i](weights, values, item_count, capacity,
                               &results[i]) &&
                    results[i].maximum_value == expected &&
                    search_knapsack_result_is_valid(
                        weights, values, item_count, capacity, &results[i]);
        }
        for (i = 0; i < 4; ++i) {
            knapsack_search_result_destroy(&results[i]);
        }
        if (!valid) {
            return 0;
        }
    }
    return 1;
}

static int ford_fulkerson_result_is_valid(
    const FordFulkersonNetwork *network,
    const FordFulkersonResult *result)
{
    int64_t *net_flow;
    size_t expected_cut_count = 0;
    size_t cut_index = 0;
    int64_t expected_cut_capacity = 0;
    size_t i;
    int valid = 1;

    if (result->vertex_count != network->vertex_count ||
        result->edge_count != network->edge_count ||
        result->source_side == NULL ||
        (network->edge_count != 0 && result->edge_flows == NULL) ||
        (result->cut_edge_count != 0 && result->cut_edges == NULL) ||
        result->source_side[result->source] == 0 ||
        result->source_side[result->sink] != 0 ||
        result->cut_capacity != result->max_flow) {
        return 0;
    }
    net_flow = (int64_t *)calloc(
        result->vertex_count, sizeof(*net_flow));
    if (net_flow == NULL) {
        return 0;
    }

    for (i = 0; i < network->edge_count; ++i) {
        const FordFulkersonEdge *edge = &network->edges[i];
        int64_t flow = result->edge_flows[i];

        if (flow < 0 || flow > edge->capacity ||
            !bellman_test_add(net_flow[edge->from], flow,
                              &net_flow[edge->from]) ||
            !bellman_test_add(net_flow[edge->to], -flow,
                              &net_flow[edge->to])) {
            valid = 0;
            break;
        }
        if (result->source_side[edge->from] != 0 &&
            result->source_side[edge->to] == 0) {
            if (cut_index >= result->cut_edge_count ||
                result->cut_edges[cut_index].from != edge->from ||
                result->cut_edges[cut_index].to != edge->to ||
                result->cut_edges[cut_index].capacity != edge->capacity ||
                !bellman_test_add(expected_cut_capacity, edge->capacity,
                                  &expected_cut_capacity)) {
                valid = 0;
                break;
            }
            ++cut_index;
            ++expected_cut_count;
        }
    }
    if (valid) {
        for (i = 0; i < result->vertex_count; ++i) {
            if (i == result->source) {
                valid = net_flow[i] == result->max_flow;
            } else if (i == result->sink) {
                valid = net_flow[i] == -result->max_flow;
            } else {
                valid = net_flow[i] == 0;
            }
            if (!valid) {
                break;
            }
        }
    }
    valid = valid && expected_cut_count == result->cut_edge_count &&
            expected_cut_capacity == result->cut_capacity;
    free(net_flow);
    return valid;
}

static int run_ford_fulkerson_case(const FordFulkersonNetwork *network,
                                   size_t source,
                                   size_t sink,
                                   int64_t expected_flow)
{
    FordFulkersonResult result = {0};
    int valid = ford_fulkerson_max_flow(
                    network, source, sink, &result) &&
                result.max_flow == expected_flow &&
                ford_fulkerson_result_is_valid(network, &result);

    ford_fulkerson_result_destroy(&result);
    return valid;
}

static int test_ford_fulkerson_flow_examples(void)
{
    const FordFulkersonEdge textbook_edges[] = {
        {0, 1, 16}, {0, 2, 13}, {1, 2, 10}, {2, 1, 4},
        {1, 3, 12}, {3, 2, 9}, {2, 4, 14}, {4, 3, 7},
        {3, 5, 20}, {4, 5, 4}
    };
    const FordFulkersonEdge disconnected_edges[] = {
        {0, 1, 5}, {2, 3, 7}
    };
    const FordFulkersonEdge zero_edges[] = {
        {0, 1, 0}, {1, 2, 9}, {0, 2, 0}
    };
    const FordFulkersonEdge multiple_paths[] = {
        {0, 1, 4}, {0, 2, 3}, {1, 3, 4}, {2, 3, 3}
    };
    const FordFulkersonEdge direct_edge[] = {{0, 1, 11}};
    const FordFulkersonEdge parallel_edges[] = {
        {0, 1, 3}, {0, 1, 4}, {1, 2, 10}
    };
    const FordFulkersonNetwork textbook = {
        6, 10, textbook_edges
    };
    const FordFulkersonNetwork disconnected = {
        4, 2, disconnected_edges
    };
    const FordFulkersonNetwork zero = {3, 3, zero_edges};
    const FordFulkersonNetwork multiple = {
        4, 4, multiple_paths
    };
    const FordFulkersonNetwork direct = {2, 1, direct_edge};
    const FordFulkersonNetwork parallel = {3, 3, parallel_edges};
    const FordFulkersonNetwork empty = {3, 0, NULL};

    return run_ford_fulkerson_case(&textbook, 0, 5, 23) &&
           run_ford_fulkerson_case(&disconnected, 0, 3, 0) &&
           run_ford_fulkerson_case(&zero, 0, 2, 0) &&
           run_ford_fulkerson_case(&multiple, 0, 3, 7) &&
           run_ford_fulkerson_case(&direct, 0, 1, 11) &&
           run_ford_fulkerson_case(&parallel, 0, 2, 7) &&
           run_ford_fulkerson_case(&empty, 0, 2, 0);
}

static int test_ford_fulkerson_reverse_residual_edges(void)
{
    const FordFulkersonEdge edges[] = {
        {0, 2, 1}, {0, 1, 1}, {1, 4, 1}, {1, 3, 1},
        {2, 3, 1}, {3, 5, 1}, {4, 5, 1}
    };
    const FordFulkersonNetwork network = {6, 7, edges};
    FordFulkersonResult result = {0};
    int valid;

    valid = ford_fulkerson_max_flow(&network, 0, 5, &result) &&
            result.max_flow == 2 &&
            result.edge_flows[0] == 1 &&
            result.edge_flows[1] == 1 &&
            result.edge_flows[2] == 1 &&
            result.edge_flows[3] == 0 &&
            result.edge_flows[4] == 1 &&
            result.edge_flows[5] == 1 &&
            result.edge_flows[6] == 1 &&
            ford_fulkerson_result_is_valid(&network, &result);
    ford_fulkerson_result_destroy(&result);
    return valid;
}

static int test_ford_fulkerson_validation_and_overflow(void)
{
    const FordFulkersonEdge valid_edges[] = {{0, 1, 5}};
    const FordFulkersonEdge negative_edges[] = {{0, 1, -1}};
    const FordFulkersonEdge invalid_edges[] = {{0, 2, 1}};
    const FordFulkersonEdge overflow_edges[] = {
        {0, 1, INT64_MAX}, {0, 2, 1}, {2, 1, 1}
    };
    const FordFulkersonNetwork valid_network = {2, 1, valid_edges};
    const FordFulkersonNetwork negative_network = {2, 1, negative_edges};
    const FordFulkersonNetwork invalid_network = {2, 1, invalid_edges};
    const FordFulkersonNetwork overflow_network = {3, 3, overflow_edges};
    const FordFulkersonNetwork null_edges = {2, 1, NULL};
    const FordFulkersonNetwork no_vertices = {0, 0, NULL};
    const FordFulkersonNetwork huge_vertices = {SIZE_MAX, 0, NULL};
    const FordFulkersonNetwork huge_edges = {2, SIZE_MAX, valid_edges};
    FordFulkersonResult result = {0};
    unsigned char *saved_source_side;
    int64_t *saved_edge_flows;
    int valid;

    if (!ford_fulkerson_max_flow(&valid_network, 0, 1, &result)) {
        return 0;
    }
    saved_source_side = result.source_side;
    saved_edge_flows = result.edge_flows;
    if (ford_fulkerson_max_flow(NULL, 0, 1, &result) ||
        ford_fulkerson_max_flow(&valid_network, 2, 1, &result) ||
        ford_fulkerson_max_flow(&valid_network, 0, 2, &result) ||
        ford_fulkerson_max_flow(&valid_network, 0, 0, &result) ||
        ford_fulkerson_max_flow(&negative_network, 0, 1, &result) ||
        ford_fulkerson_max_flow(&invalid_network, 0, 1, &result) ||
        ford_fulkerson_max_flow(&null_edges, 0, 1, &result) ||
        ford_fulkerson_max_flow(&no_vertices, 0, 0, &result) ||
        ford_fulkerson_max_flow(&huge_vertices, 0, 1, &result) ||
        ford_fulkerson_max_flow(&huge_edges, 0, 1, &result) ||
        ford_fulkerson_max_flow(&overflow_network, 0, 1, &result) ||
        ford_fulkerson_max_flow(&valid_network, 0, 1, NULL) ||
        result.source_side != saved_source_side ||
        result.edge_flows != saved_edge_flows ||
        result.max_flow != 5) {
        ford_fulkerson_result_destroy(&result);
        return 0;
    }
    valid = ford_fulkerson_result_is_valid(&valid_network, &result);
    ford_fulkerson_result_destroy(&result);
    return valid && result.source_side == NULL &&
           result.edge_flows == NULL && result.cut_edges == NULL &&
           result.max_flow == 0;
}

static int bipartite_result_is_valid(const BipartiteGraph *graph,
                                     const BipartiteMatchingResult *result)
{
    unsigned char *seen_left = NULL;
    unsigned char *seen_right = NULL;
    size_t pair_count = 0;
    size_t i;
    int valid = 1;

    if (result->left_vertex_count != graph->left_vertex_count ||
        result->right_vertex_count != graph->right_vertex_count ||
        (result->matching_size != 0 && result->pairs == NULL) ||
        (graph->left_vertex_count != 0 && result->left_to_right == NULL) ||
        (graph->right_vertex_count != 0 && result->right_to_left == NULL)) {
        return 0;
    }
    if (graph->left_vertex_count != 0) {
        seen_left = (unsigned char *)calloc(
            graph->left_vertex_count, sizeof(*seen_left));
    }
    if (graph->right_vertex_count != 0) {
        seen_right = (unsigned char *)calloc(
            graph->right_vertex_count, sizeof(*seen_right));
    }
    if ((graph->left_vertex_count != 0 && seen_left == NULL) ||
        (graph->right_vertex_count != 0 && seen_right == NULL)) {
        free(seen_left);
        free(seen_right);
        return 0;
    }

    for (i = 0; i < result->matching_size; ++i) {
        size_t left = result->pairs[i].left;
        size_t right = result->pairs[i].right;
        size_t edge_index;
        int edge_exists = 0;

        if (left >= graph->left_vertex_count ||
            right >= graph->right_vertex_count ||
            seen_left[left] != 0 || seen_right[right] != 0 ||
            result->left_to_right[left] != right ||
            result->right_to_left[right] != left ||
            (i != 0 && result->pairs[i - 1].left >= left)) {
            valid = 0;
            break;
        }
        for (edge_index = 0; edge_index < graph->edge_count; ++edge_index) {
            if (graph->edges[edge_index].left == left &&
                graph->edges[edge_index].right == right) {
                edge_exists = 1;
                break;
            }
        }
        if (!edge_exists) {
            valid = 0;
            break;
        }
        seen_left[left] = 1;
        seen_right[right] = 1;
        ++pair_count;
    }
    for (i = 0; valid && i < graph->left_vertex_count; ++i) {
        size_t right = result->left_to_right[i];

        if ((right == SIZE_MAX && seen_left[i] != 0) ||
            (right != SIZE_MAX &&
             (right >= graph->right_vertex_count ||
              result->right_to_left[right] != i || seen_left[i] == 0))) {
            valid = 0;
        }
    }
    for (i = 0; valid && i < graph->right_vertex_count; ++i) {
        size_t left = result->right_to_left[i];

        if ((left == SIZE_MAX && seen_right[i] != 0) ||
            (left != SIZE_MAX &&
             (left >= graph->left_vertex_count ||
              result->left_to_right[left] != i || seen_right[i] == 0))) {
            valid = 0;
        }
    }
    valid = valid && pair_count == result->matching_size;
    free(seen_left);
    free(seen_right);
    return valid;
}

static int run_bipartite_matching_case(const BipartiteGraph *graph,
                                       size_t expected_size)
{
    BipartiteMatchingResult result = {0};
    int valid = bipartite_maximum_matching(graph, &result) &&
                result.matching_size == expected_size &&
                bipartite_result_is_valid(graph, &result);

    bipartite_matching_result_destroy(&result);
    return valid;
}

static int test_bipartite_matching_examples(void)
{
    const BipartiteEdge perfect_edges[] = {
        {0, 0}, {0, 1}, {1, 1}, {1, 2},
        {2, 2}, {2, 3}, {3, 0}, {3, 3}
    };
    const BipartiteEdge left_heavy_edges[] = {
        {0, 0}, {1, 0}, {2, 1}, {3, 1}
    };
    const BipartiteEdge right_heavy_edges[] = {
        {0, 0}, {0, 1}, {1, 2}, {1, 3}
    };
    const BipartiteEdge reassignment_edges[] = {
        {0, 0}, {0, 1}, {1, 0}
    };
    const BipartiteEdge duplicate_edges[] = {
        {0, 1}, {0, 1}, {0, 1}, {1, 0}, {1, 0}
    };
    const BipartiteGraph perfect = {4, 4, 8, perfect_edges};
    const BipartiteGraph no_edges = {3, 4, 0, NULL};
    const BipartiteGraph left_heavy = {4, 2, 4, left_heavy_edges};
    const BipartiteGraph right_heavy = {2, 4, 4, right_heavy_edges};
    const BipartiteGraph reassignment = {2, 2, 3, reassignment_edges};
    const BipartiteGraph duplicates = {2, 2, 5, duplicate_edges};
    const BipartiteGraph empty_left = {0, 3, 0, NULL};
    const BipartiteGraph empty_right = {3, 0, 0, NULL};
    BipartiteMatchingResult result = {0};
    int valid;

    if (!run_bipartite_matching_case(&perfect, 4) ||
        !run_bipartite_matching_case(&no_edges, 0) ||
        !run_bipartite_matching_case(&left_heavy, 2) ||
        !run_bipartite_matching_case(&right_heavy, 2) ||
        !run_bipartite_matching_case(&duplicates, 2) ||
        !run_bipartite_matching_case(&empty_left, 0) ||
        !run_bipartite_matching_case(&empty_right, 0) ||
        !bipartite_maximum_matching(&reassignment, &result)) {
        bipartite_matching_result_destroy(&result);
        return 0;
    }
    valid = result.matching_size == 2 &&
            result.left_to_right[0] == 1 &&
            result.left_to_right[1] == 0 &&
            bipartite_result_is_valid(&reassignment, &result);
    bipartite_matching_result_destroy(&result);
    return valid;
}

static int test_bipartite_matching_validation(void)
{
    const BipartiteEdge valid_edges[] = {{0, 0}, {1, 1}};
    const BipartiteEdge invalid_left_edges[] = {{2, 0}};
    const BipartiteEdge invalid_right_edges[] = {{0, 2}};
    const BipartiteGraph valid_graph = {2, 2, 2, valid_edges};
    const BipartiteGraph invalid_left = {2, 2, 1, invalid_left_edges};
    const BipartiteGraph invalid_right = {2, 2, 1, invalid_right_edges};
    const BipartiteGraph null_edges = {2, 2, 1, NULL};
    const BipartiteGraph huge_left = {SIZE_MAX, 0, 0, NULL};
    const BipartiteGraph huge_right = {0, SIZE_MAX, 0, NULL};
    const BipartiteGraph huge_edges = {2, 2, SIZE_MAX, valid_edges};
    BipartiteMatchingResult result = {0};
    size_t *saved_left_map;
    size_t *saved_right_map;
    int valid;

    if (!bipartite_maximum_matching(&valid_graph, &result)) {
        return 0;
    }
    saved_left_map = result.left_to_right;
    saved_right_map = result.right_to_left;
    if (bipartite_maximum_matching(NULL, &result) ||
        bipartite_maximum_matching(&invalid_left, &result) ||
        bipartite_maximum_matching(&invalid_right, &result) ||
        bipartite_maximum_matching(&null_edges, &result) ||
        bipartite_maximum_matching(&huge_left, &result) ||
        bipartite_maximum_matching(&huge_right, &result) ||
        bipartite_maximum_matching(&huge_edges, &result) ||
        bipartite_maximum_matching(&valid_graph, NULL) ||
        result.left_to_right != saved_left_map ||
        result.right_to_left != saved_right_map ||
        result.matching_size != 2) {
        bipartite_matching_result_destroy(&result);
        return 0;
    }
    valid = bipartite_result_is_valid(&valid_graph, &result);
    bipartite_matching_result_destroy(&result);
    return valid && result.pairs == NULL &&
           result.left_to_right == NULL && result.right_to_left == NULL &&
           result.matching_size == 0;
}

static void brute_force_matching_search(const BipartiteGraph *graph,
                                        size_t left,
                                        uint64_t used_right,
                                        size_t current_size,
                                        size_t *maximum_size)
{
    size_t edge_index;

    if (left == graph->left_vertex_count) {
        if (current_size > *maximum_size) {
            *maximum_size = current_size;
        }
        return;
    }

    brute_force_matching_search(graph, left + 1, used_right,
                                current_size, maximum_size);
    for (edge_index = 0; edge_index < graph->edge_count; ++edge_index) {
        size_t right = graph->edges[edge_index].right;
        uint64_t bit;

        if (graph->edges[edge_index].left != left || right >= 64) {
            continue;
        }
        bit = UINT64_C(1) << right;
        if ((used_right & bit) == 0) {
            brute_force_matching_search(graph, left + 1,
                                        used_right | bit,
                                        current_size + 1, maximum_size);
        }
    }
}

static int test_bipartite_matching_brute_force_agreement(void)
{
    uint32_t random_state = UINT32_C(0xb17a4e);
    size_t test_case;

    for (test_case = 0; test_case < 64; ++test_case) {
        BipartiteEdge edges[72];
        size_t left_count =
            (size_t)(next_knapsack_test_random(&random_state) % 7);
        size_t right_count =
            (size_t)(next_knapsack_test_random(&random_state) % 7);
        size_t edge_count = 0;
        size_t left;
        size_t right;
        size_t expected = 0;
        BipartiteGraph graph;
        BipartiteMatchingResult result = {0};
        size_t i;
        int valid;

        for (left = 0; left < left_count; ++left) {
            for (right = 0; right < right_count; ++right) {
                if (next_knapsack_test_random(&random_state) % 3 == 0) {
                    edges[edge_count].left = left;
                    edges[edge_count].right = right;
                    ++edge_count;
                    if (next_knapsack_test_random(&random_state) % 4 == 0) {
                        edges[edge_count] = edges[edge_count - 1];
                        ++edge_count;
                    }
                }
            }
        }
        graph.left_vertex_count = left_count;
        graph.right_vertex_count = right_count;
        graph.edge_count = edge_count;
        graph.edges = edge_count == 0 ? NULL : edges;
        brute_force_matching_search(&graph, 0, 0, 0, &expected);
        valid = bipartite_maximum_matching(&graph, &result) &&
                result.matching_size == expected &&
                bipartite_result_is_valid(&graph, &result);
        for (i = 0; valid && i < result.matching_size; ++i) {
            size_t edge_index;
            int exists = 0;

            for (edge_index = 0; edge_index < edge_count; ++edge_index) {
                if (edges[edge_index].left == result.pairs[i].left &&
                    edges[edge_index].right == result.pairs[i].right) {
                    exists = 1;
                    break;
                }
            }
            valid = exists;
        }
        bipartite_matching_result_destroy(&result);
        if (!valid) {
            return 0;
        }
    }
    return 1;
}

static int test_fast_exponentiation(void)
{
    double result;
    DcaOperationCounts operations;

    if (!dca_fast_power(2.0, 10, &result, &operations) ||
        result != 1024.0 || operations.multiplications == 0) {
        return 0;
    }
    if (!dca_fast_power(0.0, 0, &result, NULL) || result != 1.0 ||
        !dca_fast_power(2.0, -3, &result, NULL) || result != 0.125 ||
        !dca_fast_power(-2.0, 3, &result, NULL) || result != -8.0 ||
        !dca_fast_power(1.0, INT_MIN, &result, NULL) || result != 1.0) {
        return 0;
    }
    if (dca_fast_power(0.0, -1, &result, NULL) ||
        dca_fast_power(DBL_MAX, 2, &result, NULL) ||
        dca_fast_power(INFINITY, 2, &result, NULL) ||
        dca_fast_power(2.0, 3, NULL, NULL)) {
        return 0;
    }

    return 1;
}

static int test_celebrity_problem(void)
{
    const unsigned char known_celebrity[] = {
        0, 1, 1, 0,
        0, 0, 1, 0,
        0, 0, 0, 0,
        0, 1, 1, 0
    };
    const unsigned char no_celebrity[] = {
        0, 0, 0,
        0, 0, 0,
        0, 0, 0
    };
    const unsigned char one_person[] = {1};
    unsigned char larger[64 * 64] = {0};
    DcaOperationCounts operations;
    size_t person;

    if (dca_find_celebrity(known_celebrity, 4, &operations) != 2 ||
        operations.comparisons == 0 ||
        dca_find_celebrity(no_celebrity, 3, NULL) != -1 ||
        dca_find_celebrity(one_person, 1, NULL) != 0 ||
        dca_find_celebrity(NULL, 1, NULL) != -1 ||
        dca_find_celebrity(one_person, 0, NULL) != -1 ||
        dca_find_celebrity(one_person, SIZE_MAX, NULL) != -1) {
        return 0;
    }

    for (person = 0; person < 64; ++person) {
        if (person != 17) {
            larger[person * 64 + 17] = 1;
        }
    }
    return dca_find_celebrity(larger, 64, &operations) == 17 &&
           operations.comparisons <= 3 * 64 - 3;
}

static int test_classical_matrix_multiplication(void)
{
    const double expected[] = {58.0, 64.0, 139.0, 154.0};
    DcaMatrix left = {0, 0, NULL};
    DcaMatrix right = {0, 0, NULL};
    DcaMatrix result = {0, 0, NULL};
    DcaMatrix zero_left = {0, 0, NULL};
    DcaMatrix zero_right = {0, 0, NULL};
    DcaMatrix zero_result = {0, 0, NULL};
    DcaMatrix overflow_left = {0, 0, NULL};
    DcaMatrix overflow_right = {0, 0, NULL};
    DcaMatrix overflow_result = {0, 0, NULL};
    DcaMatrix too_large = {0, 0, NULL};
    DcaOperationCounts operations;
    int success = 0;

    if (!dca_matrix_create(&left, 2, 3) ||
        !dca_matrix_create(&right, 3, 2)) {
        goto cleanup;
    }
    left.values[0] = 1.0;
    left.values[1] = 2.0;
    left.values[2] = 3.0;
    left.values[3] = 4.0;
    left.values[4] = 5.0;
    left.values[5] = 6.0;
    right.values[0] = 7.0;
    right.values[1] = 8.0;
    right.values[2] = 9.0;
    right.values[3] = 10.0;
    right.values[4] = 11.0;
    right.values[5] = 12.0;

    if (!dca_matrix_multiply_classical(&left, &right, &result, &operations) ||
        !matrix_matches_values(&result, 2, 2, expected) ||
        operations.multiplications != 12 || operations.additions != 8) {
        goto cleanup;
    }

    dca_matrix_destroy(&result);
    if (!dca_matrix_create(&zero_left, 2, 0) ||
        !dca_matrix_create(&zero_right, 0, 3) ||
        !dca_matrix_multiply_classical(&zero_left, &zero_right,
                                       &zero_result, NULL)) {
        goto cleanup;
    }
    {
        const double zero_values[] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
        if (!matrix_matches_values(&zero_result, 2, 3, zero_values)) {
            goto cleanup;
        }
    }

    if (dca_matrix_multiply_classical(NULL, &right, &result, NULL) ||
        dca_matrix_multiply_classical(&left, &zero_right, &result, NULL) ||
        dca_matrix_create(&too_large, SIZE_MAX, 2) ||
        dca_matrix_create(NULL, 1, 1)) {
        goto cleanup;
    }
    if (!dca_matrix_create(&overflow_left, 1, 1) ||
        !dca_matrix_create(&overflow_right, 1, 1)) {
        goto cleanup;
    }
    overflow_left.values[0] = DBL_MAX;
    overflow_right.values[0] = 2.0;
    if (dca_matrix_multiply_classical(&overflow_left, &overflow_right,
                                      &overflow_result, NULL)) {
        goto cleanup;
    }
    success = 1;

cleanup:
    dca_matrix_destroy(&left);
    dca_matrix_destroy(&right);
    dca_matrix_destroy(&result);
    dca_matrix_destroy(&zero_left);
    dca_matrix_destroy(&zero_right);
    dca_matrix_destroy(&zero_result);
    dca_matrix_destroy(&overflow_left);
    dca_matrix_destroy(&overflow_right);
    dca_matrix_destroy(&overflow_result);
    dca_matrix_destroy(&too_large);
    return success;
}

static int test_strassen_matrix_multiplication(void)
{
    DcaMatrix left = {0, 0, NULL};
    DcaMatrix right = {0, 0, NULL};
    DcaMatrix classical = {0, 0, NULL};
    DcaMatrix strassen = {0, 0, NULL};
    DcaMatrix larger_left = {0, 0, NULL};
    DcaMatrix larger_right = {0, 0, NULL};
    DcaMatrix larger_classical = {0, 0, NULL};
    DcaMatrix larger_strassen = {0, 0, NULL};
    DcaMatrix invalid_left = {0, 0, NULL};
    DcaMatrix invalid_right = {0, 0, NULL};
    DcaMatrix overflow_left = {0, 0, NULL};
    DcaMatrix overflow_right = {0, 0, NULL};
    DcaMatrix overflow_result = {0, 0, NULL};
    DcaOperationCounts operations;
    size_t i;
    int success = 0;

    if (!dca_matrix_create(&left, 3, 3) ||
        !dca_matrix_create(&right, 3, 3)) {
        goto cleanup;
    }
    for (i = 0; i < 9; ++i) {
        left.values[i] = (double)((int)i - 4);
        right.values[i] = (double)((int)(8 - i) - 2);
    }
    if (!dca_matrix_multiply_classical(&left, &right, &classical, NULL) ||
        !dca_matrix_multiply_strassen(&left, &right, &strassen, &operations) ||
        !matrices_are_close(&classical, &strassen)) {
        goto cleanup;
    }

    if (!dca_matrix_create(&larger_left, 65, 65) ||
        !dca_matrix_create(&larger_right, 65, 65)) {
        goto cleanup;
    }
    for (i = 0; i < 65 * 65; ++i) {
        larger_left.values[i] = (double)((int)(i * 7 % 11) - 5);
        larger_right.values[i] = (double)((int)(i * 3 % 13) - 6);
    }
    if (!dca_matrix_multiply_classical(&larger_left, &larger_right,
                                       &larger_classical, NULL) ||
        !dca_matrix_multiply_strassen(&larger_left, &larger_right,
                                      &larger_strassen, &operations) ||
        !matrices_are_close(&larger_classical, &larger_strassen) ||
        operations.multiplications == 0) {
        goto cleanup;
    }

    if (!dca_matrix_create(&invalid_left, 2, 3) ||
        !dca_matrix_create(&invalid_right, 2, 2) ||
        dca_matrix_multiply_strassen(NULL, &invalid_right,
                                     &overflow_result, NULL) ||
        dca_matrix_multiply_strassen(&invalid_left, &invalid_right,
                                     &overflow_result, NULL) ||
        !dca_matrix_create(&overflow_left, 1, 1) ||
        !dca_matrix_create(&overflow_right, 1, 1)) {
        goto cleanup;
    }
    overflow_left.values[0] = DBL_MAX;
    overflow_right.values[0] = 2.0;
    if (dca_matrix_multiply_strassen(&overflow_left, &overflow_right,
                                     &overflow_result, NULL)) {
        goto cleanup;
    }
    success = 1;

cleanup:
    dca_matrix_destroy(&left);
    dca_matrix_destroy(&right);
    dca_matrix_destroy(&classical);
    dca_matrix_destroy(&strassen);
    dca_matrix_destroy(&larger_left);
    dca_matrix_destroy(&larger_right);
    dca_matrix_destroy(&larger_classical);
    dca_matrix_destroy(&larger_strassen);
    dca_matrix_destroy(&invalid_left);
    dca_matrix_destroy(&invalid_right);
    dca_matrix_destroy(&overflow_left);
    dca_matrix_destroy(&overflow_right);
    dca_matrix_destroy(&overflow_result);
    return success;
}

static int sequence_records_equal(const DnaSequence *left,
                                  const DnaSequence *right)
{
    return left->length == right->length &&
           strcmp(left->identifier, right->identifier) == 0;
}

static int records_are_sorted(const DnaSequence *records,
                              size_t count,
                              SequenceSortCriterion criterion)
{
    size_t i;

    for (i = 1; i < count; ++i) {
        int out_of_order;

        if (criterion == SEQUENCE_SORT_BY_LENGTH) {
            out_of_order = records[i - 1].length > records[i].length;
        } else {
            out_of_order =
                strcmp(records[i - 1].identifier, records[i].identifier) > 0;
        }
        if (out_of_order) {
            return 0;
        }
    }

    return 1;
}

static int records_have_same_contents(const DnaSequence *actual,
                                      const DnaSequence *expected,
                                      size_t count)
{
    int matched[48] = {0};
    size_t i;

    for (i = 0; i < count; ++i) {
        size_t j;

        for (j = 0; j < count; ++j) {
            if (!matched[j] &&
                sequence_records_equal(&actual[i], &expected[j])) {
                matched[j] = 1;
                break;
            }
        }
        if (j == count) {
            return 0;
        }
    }

    return 1;
}

static int run_sort_suite(SortFunction sort,
                          SequenceSortCriterion criterion)
{
    DnaSequence source[] = {
        {"alpha", "A", 1}, {"bravo", "AA", 2},
        {"charlie", "AA", 2}, {"delta", "AAA", 3},
        {"echo", "AAA", 3}, {"foxtrot", "AAAA", 4},
        {"golf", "AAAAA", 5}, {"hotel", "AAAAAA", 6}
    };
    const size_t sorted[] = {0, 1, 2, 3, 4, 5, 6, 7};
    const size_t reversed[] = {7, 6, 5, 4, 3, 2, 1, 0};
    const size_t duplicates[] = {3, 2, 4, 2, 0, 4, 7, 3};
    const size_t randomized[] = {4, 0, 7, 2, 5, 1, 6, 3};
    const size_t *orders[] = {NULL, sorted, sorted, reversed,
                              duplicates, randomized};
    const size_t sizes[] = {0, 1, 8, 8, 8, 8};
    size_t case_index;

    for (case_index = 0; case_index < sizeof(sizes) / sizeof(sizes[0]);
         ++case_index) {
        DnaSequence input[8];
        DnaSequence original[8];
        SequenceAlgorithmStats stats;
        size_t i;

        for (i = 0; i < sizes[case_index]; ++i) {
            input[i] = source[orders[case_index][i]];
            original[i] = input[i];
        }

        if (!sort(sizes[case_index] == 0 ? NULL : input,
                  sizes[case_index], criterion, &stats) ||
            !records_are_sorted(input, sizes[case_index], criterion) ||
            !records_have_same_contents(input, original, sizes[case_index])) {
            return 0;
        }
        if (sizes[case_index] > 1 && stats.comparisons == 0) {
            return 0;
        }
    }

    {
        DnaSequence random_records[48];
        DnaSequence original[48];
        char identifiers[48][20];
        unsigned int random_state = 0xC17U;
        size_t i;

        for (i = 0; i < 48; ++i) {
            (void)snprintf(identifiers[i], sizeof(identifiers[i]),
                           "record-%02lu", (unsigned long)i);
            random_state = random_state * 1664525U + 1013904223U;
            random_records[i].identifier = identifiers[i];
            random_records[i].bases = "A";
            random_records[i].length = random_state % 13U;
        }
        for (i = 47; i > 0; --i) {
            size_t other;
            DnaSequence temporary;

            random_state = random_state * 1664525U + 1013904223U;
            other = random_state % (i + 1);
            temporary = random_records[i];
            random_records[i] = random_records[other];
            random_records[other] = temporary;
        }
        for (i = 0; i < 48; ++i) {
            original[i] = random_records[i];
        }

        {
            SequenceAlgorithmStats stats;

            if (!sort(random_records, 48, criterion, &stats) ||
                !records_are_sorted(random_records, 48, criterion) ||
                !records_have_same_contents(random_records, original, 48) ||
                stats.comparisons == 0) {
                return 0;
            }
        }
    }

    return 1;
}

static int run_search_suite(SearchFunction search,
                            SequenceSortCriterion criterion)
{
    DnaSequence records[] = {
        {"alpha", "A", 1}, {"bravo", "AA", 2},
        {"charlie", "AA", 2}, {"delta", "AAA", 3},
        {"echo", "AAA", 3}, {"foxtrot", "AAAA", 4},
        {"golf", "AAAAA", 5}, {"hotel", "AAAAAA", 6}
    };
    DnaSequence first = {"alpha", "A", 1};
    DnaSequence last = {"hotel", "AAAAAA", 6};
    DnaSequence duplicate_length = {NULL, NULL, 3};
    DnaSequence missing_identifier = {"zulu", "A", 1};
    DnaSequence missing_length = {NULL, NULL, 7};
    DnaSequence single = {"only", "AA", 2};
    size_t index;
    SequenceAlgorithmStats stats;
    const DnaSequence *missing =
        criterion == SEQUENCE_SORT_BY_LENGTH
            ? &missing_length
            : &missing_identifier;

    if (search(records, 8, &first, criterion, &index, &stats) !=
            SEQUENCE_SEARCH_FOUND ||
        index != 0 || stats.comparisons == 0) {
        return 0;
    }
    if (search(records, 8, &last, criterion, &index, &stats) !=
            SEQUENCE_SEARCH_FOUND ||
        index != 7) {
        return 0;
    }
    if (search(records, 8, missing, criterion, &index, &stats) !=
            SEQUENCE_SEARCH_NOT_FOUND ||
        index != SIZE_MAX) {
        return 0;
    }
    if (criterion == SEQUENCE_SORT_BY_LENGTH) {
        if (search(records, 8, &duplicate_length, criterion, &index,
                   &stats) != SEQUENCE_SEARCH_FOUND ||
            records[index].length != duplicate_length.length) {
            return 0;
        }
    }
    if (search(&single, 1, &single, criterion, &index, &stats) !=
            SEQUENCE_SEARCH_FOUND ||
        index != 0) {
        return 0;
    }
    if (search(NULL, 0, &first, criterion, &index, &stats) !=
            SEQUENCE_SEARCH_NOT_FOUND ||
        index != SIZE_MAX) {
        return 0;
    }

    return 1;
}

static int test_merge_sort_by_length(void)
{
    return run_sort_suite(sequence_merge_sort, SEQUENCE_SORT_BY_LENGTH);
}

static int test_merge_sort_by_identifier(void)
{
    return run_sort_suite(sequence_merge_sort, SEQUENCE_SORT_BY_IDENTIFIER);
}

static int test_quick_sort_by_length(void)
{
    return run_sort_suite(sequence_quick_sort, SEQUENCE_SORT_BY_LENGTH);
}

static int test_quick_sort_by_identifier(void)
{
    return run_sort_suite(sequence_quick_sort, SEQUENCE_SORT_BY_IDENTIFIER);
}

static int test_iterative_binary_search_by_length(void)
{
    return run_search_suite(sequence_binary_search_iterative,
                            SEQUENCE_SORT_BY_LENGTH);
}

static int test_iterative_binary_search_by_identifier(void)
{
    return run_search_suite(sequence_binary_search_iterative,
                            SEQUENCE_SORT_BY_IDENTIFIER);
}

static int test_recursive_binary_search_by_length(void)
{
    return run_search_suite(sequence_binary_search_recursive,
                            SEQUENCE_SORT_BY_LENGTH);
}

static int test_recursive_binary_search_by_identifier(void)
{
    return run_search_suite(sequence_binary_search_recursive,
                            SEQUENCE_SORT_BY_IDENTIFIER);
}

static int test_valid_sequences(void)
{
    return dna_sequence_is_valid("ACGTacgt") &&
           !dna_sequence_is_valid("ACGTX") &&
           !dna_sequence_is_valid("") &&
           !dna_sequence_is_valid(NULL);
}

static int test_uppercase_sequence(void)
{
    return dna_sequence_is_valid("ACGT");
}

static int test_lowercase_sequence(void)
{
    return dna_sequence_is_valid("acgt");
}

static int test_single_nucleotide(void)
{
    return dna_sequence_is_valid("G");
}

static int test_invalid_characters(void)
{
    return !dna_sequence_is_valid("ACGTN") &&
           !dna_sequence_is_valid("ACG-T");
}

static int test_empty_sequence(void)
{
    return !dna_sequence_is_valid("");
}

static int test_null_sequence(void)
{
    return !dna_sequence_is_valid(NULL);
}

static int test_create_and_destroy(void)
{
    DnaSequence sequence = {NULL, NULL, 0};

    if (!dna_sequence_create(&sequence, "sample-1", "ACGT")) {
        return 0;
    }

    if (strcmp(sequence.identifier, "sample-1") != 0 ||
        strcmp(sequence.bases, "ACGT") != 0 ||
        sequence.length != 4) {
        dna_sequence_destroy(&sequence);
        return 0;
    }

    dna_sequence_destroy(&sequence);
    dna_sequence_destroy(&sequence);

    return sequence.identifier == NULL &&
           sequence.bases == NULL &&
           sequence.length == 0;
}

static int test_create_into_occupied_structure(void)
{
    DnaSequence sequence = {NULL, NULL, 0};
    char *identifier;
    char *bases;
    int result;
    int unchanged;

    if (!dna_sequence_create(&sequence, "original", "ACGT"))
    {
        return 0;
    }

    identifier = sequence.identifier;
    bases = sequence.bases;
    result = dna_sequence_create(&sequence, "replacement", "TGCA");
    unchanged = sequence.identifier == identifier &&
                sequence.bases == bases &&
                sequence.length == 4 &&
                strcmp(sequence.identifier, "original") == 0 &&
                strcmp(sequence.bases, "ACGT") == 0;

    dna_sequence_destroy(&sequence);
    return !result && unchanged;
}

static int test_invalid_identifier(void)
{
    DnaSequence sequence = {NULL, NULL, 0};
    int empty_rejected;
    int null_rejected;

    empty_rejected = !dna_sequence_create(&sequence, "", "ACGT") &&
                     sequence.identifier == NULL &&
                     sequence.bases == NULL &&
                     sequence.length == 0;
    null_rejected = !dna_sequence_create(&sequence, NULL, "ACGT") &&
                    sequence.identifier == NULL &&
                    sequence.bases == NULL &&
                    sequence.length == 0;

    return empty_rejected && null_rejected;
}

static int test_invalid_sequence_creation(void)
{
    DnaSequence sequence = {NULL, NULL, 0};

    return !dna_sequence_create(&sequence, "sample", "ACGX") &&
           sequence.identifier == NULL &&
           sequence.bases == NULL &&
           sequence.length == 0;
}

static int write_fixture(const char *path,
                         const char *contents,
                         size_t contents_length)
{
    FILE *file = fopen(path, "wb");
    int success;

    if (file == NULL) {
        return 0;
    }

    success = fwrite(contents, 1, contents_length, file) == contents_length;
    if (fclose(file) != 0) {
        success = 0;
    }
    return success;
}

static int load_fixture(const char *path,
                        const char *contents,
                        size_t contents_length,
                        DnaSequence **sequences,
                        size_t *count)
{
    int success;

    *sequences = NULL;
    *count = 0;
    if (!write_fixture(path, contents, contents_length)) {
        remove(path);
        return 0;
    }

    success = fasta_read_file(path, sequences, count);
    remove(path);
    return success;
}

static int test_fasta_single_record(void)
{
    const char contents[] = ">sample\nACGT\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_single.fasta", contents,
                               sizeof(contents) - 1, &sequences, &count);
    int valid = success && count == 1 &&
                strcmp(sequences[0].identifier, "sample") == 0 &&
                strcmp(sequences[0].bases, "ACGT") == 0 &&
                sequences[0].length == 4;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_multiple_records(void)
{
    const char contents[] = ">first\nAC\n>second description\nGT\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_multiple.fasta", contents,
                               sizeof(contents) - 1, &sequences, &count);
    int valid = success && count == 2 &&
                strcmp(sequences[0].identifier, "first") == 0 &&
                strcmp(sequences[0].bases, "AC") == 0 &&
                strcmp(sequences[1].identifier, "second description") == 0 &&
                strcmp(sequences[1].bases, "GT") == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_multiline_sequence(void)
{
    const char contents[] = ">multiline\nAC\nGT\nac\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_multiline.fasta", contents,
                               sizeof(contents) - 1, &sequences, &count);
    int valid = success && count == 1 &&
                strcmp(sequences[0].bases, "ACGTac") == 0 &&
                sequences[0].length == 6;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_invalid_dna(void)
{
    const char contents[] = ">invalid\nACNX\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_invalid_dna.fasta", contents,
                               sizeof(contents) - 1, &sequences, &count);
    int valid = !success && sequences == NULL && count == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_missing_header(void)
{
    const char contents[] = "ACGT\n>valid\nTGCA\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_missing_header.fasta", contents,
                               sizeof(contents) - 1, &sequences, &count);
    int valid = !success && sequences == NULL && count == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_empty_identifier(void)
{
    const char contents[] = ">\nACGT\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_empty_identifier.fasta", contents,
                               sizeof(contents) - 1, &sequences, &count);
    int valid = !success && sequences == NULL && count == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_whitespace_identifier(void)
{
    const char contents[] = "> \t \nACGT\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_whitespace_identifier.fasta",
                               contents, sizeof(contents) - 1,
                               &sequences, &count);
    int valid = !success && sequences == NULL && count == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_embedded_null(void)
{
    const char contents[] = {
        '>', 'n', 'u', 'l', 'l', '\n', 'A', 'C', '\0', 'G', 'T', '\n'
    };
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_embedded_null.fasta", contents,
                               sizeof(contents), &sequences, &count);
    int valid = !success && sequences == NULL && count == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_empty_sequence(void)
{
    const char contents[] = ">empty\n>next\nAC\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_empty_sequence.fasta", contents,
                               sizeof(contents) - 1, &sequences, &count);
    int valid = !success && sequences == NULL && count == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_empty_file(void)
{
    const char contents[] = "";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_empty_file.fasta", contents, 0,
                               &sequences, &count);
    int valid = !success && sequences == NULL && count == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_missing_file(void)
{
    DnaSequence *sequences = NULL;
    size_t count = 27;

    remove("genome_test_missing_file.fasta");
    return !fasta_read_file("genome_test_missing_file.fasta",
                            &sequences, &count) &&
           sequences == NULL && count == 0;
}

static int test_fasta_null_output_parameters(void)
{
    DnaSequence dummy = {NULL, NULL, 0};
    DnaSequence *sequences = &dummy;
    size_t count = 42;
    int null_sequences_rejected;
    int null_count_rejected;

    null_sequences_rejected =
        !fasta_read_file("", NULL, &count) && count == 0;
    null_count_rejected =
        !fasta_read_file("", &sequences, NULL) && sequences == NULL;
    return null_sequences_rejected && null_count_rejected;
}

static int test_fasta_crlf(void)
{
    const char contents[] = ">windows\r\nAC\r\nGT\r\n";
    DnaSequence *sequences = NULL;
    size_t count = 0;
    int success = load_fixture("genome_test_crlf.fasta", contents,
                               sizeof(contents) - 1, &sequences, &count);
    int valid = success && count == 1 &&
                strcmp(sequences[0].identifier, "windows") == 0 &&
                strcmp(sequences[0].bases, "ACGT") == 0;

    fasta_sequences_destroy(sequences, count);
    return valid;
}

static int test_fasta_long_sequence(void)
{
    const size_t sequence_length = 16384;
    const size_t line_length = 1024;
    const size_t line_count = sequence_length / line_length;
    const char header[] = ">long\n";
    size_t contents_size = sizeof(header) - 1 + sequence_length +
                           line_count + 1;
    char *contents = (char *)malloc(contents_size);
    DnaSequence *sequences = NULL;
    size_t count = 0;
    size_t offset = sizeof(header) - 1;
    size_t i;
    int success;
    int valid;

    if (contents == NULL) {
        return 0;
    }

    memcpy(contents, header, sizeof(header) - 1);
    for (i = 0; i < line_count; ++i)
    {
        memset(contents + offset, 'A', line_length);
        offset += line_length;
        contents[offset++] = '\n';
    }
    contents[contents_size - 1] = '\0';

    success = load_fixture("genome_test_long.fasta", contents,
                           contents_size - 1, &sequences, &count);
    valid = success && count == 1 &&
            sequences[0].length == sequence_length &&
            strlen(sequences[0].bases) == sequence_length;

    fasta_sequences_destroy(sequences, count);
    free(contents);
    return valid;
}

static int test_greedy_activity_selection(void)
{
    Activity activities[] = {
        {1, 4, 1}, {3, 5, 2}, {0, 6, 3},
        {5, 7, 4}, {8, 9, 5}, {5, 9, 6}
    };
    Activity selected[6];
    size_t count = 0;

    return greedy_activity_selection(activities, 6, selected, 6, &count) &&
           count == 3 && selected[0].identifier == 1 &&
           selected[1].identifier == 4 && selected[2].identifier == 5 &&
           greedy_activity_selection(NULL, 0, NULL, 0, &count) &&
           count == 0 &&
           !greedy_activity_selection(
               (Activity[]){{4, 4, 1}}, 1, NULL, 0, &count);
}

static int test_greedy_fractional_knapsack(void)
{
    KnapsackItem items[] = {
        {10.0, 60.0, 1}, {20.0, 100.0, 2},
        {30.0, 120.0, 3}, {0.0, 25.0, 4}
    };
    double fractions[4] = {0};
    double total = 0.0;

    return greedy_fractional_knapsack(items, 4, 50.0,
                                      fractions, &total) &&
           fabs(total - 265.0) < 1e-9 &&
           fractions[0] == 10.0 && fractions[1] == 20.0 &&
           fractions[2] == 20.0 &&
           !greedy_fractional_knapsack(items, 4, -1.0, NULL, NULL) &&
           !greedy_fractional_knapsack(
               (KnapsackItem[]){{-1.0, 10.0, 1}}, 1, 1.0, NULL, NULL);
}

static int test_greedy_huffman(void)
{
    const unsigned char symbols[] = {'A', 'B', 'C', 'D'};
    const double frequencies[] = {5.0, 2.0, 1.0, 1.0};
    const unsigned char input[] = {'A', 'B', 'C', 'D', 'A'};
    unsigned char *bits = NULL, *decoded = NULL;
    size_t bit_count = 0, decoded_count = 0;
    HuffmanNode *root = greedy_huffman_build(symbols, frequencies, 4);
    HuffmanNode *single = greedy_huffman_build(
        (const unsigned char *)"A", (const double[]) {1.0}, 1);
    unsigned char *single_bits = NULL, *single_decoded = NULL;
    size_t single_bit_count = 0, single_decoded_count = 0;
    int valid;

    if (root == NULL) return 0;
    valid = greedy_huffman_encode(input, 5, root, &bits, &bit_count) &&
            greedy_huffman_decode(bits, bit_count, root,
                                  &decoded, &decoded_count) &&
            decoded_count == 5 &&
            memcmp(decoded, input, 5) == 0;
    free(bits);
    free(decoded);
    bits = NULL;
    decoded = NULL;
    if (valid) {
        valid = !greedy_huffman_decode(
            (const unsigned char[]){2}, 1, root, &decoded, &decoded_count);
    }
    greedy_huffman_destroy(root);

    if (!valid || single == NULL) {
        greedy_huffman_destroy(single);
        return 0;
    }
    valid = greedy_huffman_encode(
                (const unsigned char *)"A", 1, single,
                &single_bits, &single_bit_count) &&
            single_bit_count == 1 &&
            greedy_huffman_decode(single_bits, single_bit_count, single,
                                  &single_decoded, &single_decoded_count) &&
            single_decoded_count == 1 &&
            single_decoded[0] == 'A';
    free(single_bits);
    free(single_decoded);
    single_decoded = NULL;
    single_decoded_count = 0;
    if (valid) {
        valid = !greedy_huffman_decode(
                    (const unsigned char[]){0}, SIZE_MAX, single,
                    &single_decoded, &single_decoded_count) &&
                single_decoded == NULL && single_decoded_count == 0;
    }
    greedy_huffman_destroy(single);
    return valid;
}

static int test_greedy_huffman_validation(void)
{
    unsigned char symbols[257];
    size_t i;
    const unsigned char duplicate[] = {'A', 'A'};
    const double normal[] = {1.0, 2.0};
    const double negative[] = {1.0, -1.0};
    const double not_a_number[] = {1.0, NAN};
    const double infinite[] = {1.0, INFINITY};
    const double overflowing[] = {DBL_MAX, DBL_MAX};

    for (i = 0; i < sizeof(symbols); ++i) {
        symbols[i] = (unsigned char)i;
    }
    return greedy_huffman_build(duplicate, normal, 2) == NULL &&
           greedy_huffman_build(symbols, normal, 257) == NULL &&
           greedy_huffman_build((const unsigned char *)"AB", negative, 2) ==
               NULL &&
           greedy_huffman_build((const unsigned char *)"AB",
                                not_a_number, 2) == NULL &&
           greedy_huffman_build((const unsigned char *)"AB", infinite, 2) ==
               NULL &&
           greedy_huffman_build((const unsigned char *)"AB", overflowing,
                                2) == NULL;
}

static int test_greedy_huffman_edge_cases(void)
{
    const unsigned char symbols[] = {'A', 'B'};
    const double frequencies[] = {1.0, 1.0};
    const unsigned char partial_invalid_bits[] = {0, 2};
    const unsigned char unknown[] = {'C'};
    HuffmanNode *root = greedy_huffman_build(symbols, frequencies, 2);
    unsigned char *bits = NULL;
    unsigned char *decoded = NULL;
    size_t bit_count = 0, decoded_count = 0;
    int valid;

    if (root == NULL) return 0;
    valid = !greedy_huffman_decode(partial_invalid_bits, 2, root,
                                   &decoded, &decoded_count) &&
            decoded == NULL && decoded_count == 0;
    if (valid) {
        valid = !greedy_huffman_encode(unknown, 1, root, &bits, &bit_count) &&
                bits == NULL && bit_count == 0;
    }
    if (valid) {
        valid = greedy_huffman_encode(NULL, 0, root, &bits, &bit_count) &&
                bits != NULL && bit_count == 0 &&
                greedy_huffman_decode(bits, bit_count, root,
                                      &decoded, &decoded_count) &&
                decoded != NULL && decoded_count == 0;
    }
    free(bits);
    free(decoded);
    greedy_huffman_destroy(root);
    return valid;
}

static int test_greedy_dijkstra(void)
{
    double adjacency[16] = {
        0.0, 1.0, DBL_MAX, 4.0,
        1.0, 0.0, 2.0, DBL_MAX,
        DBL_MAX, 2.0, 0.0, 1.0,
        4.0, DBL_MAX, 1.0, 0.0
    };
    WeightedGraph graph = {4, adjacency};
    ShortestPathResult result = {0, NULL, NULL};
    int valid;

    valid = greedy_dijkstra_matrix(&graph, 0, &result) &&
            result.distances[0] == 0.0 &&
            result.distances[1] == 1.0 &&
            result.distances[2] == 3.0 &&
            result.distances[3] == 4.0 &&
            result.predecessors[2] == 1 &&
            !greedy_dijkstra_matrix(&graph, 4, &result);
    greedy_shortest_path_destroy(&result);
    return valid;
}

static int test_dijkstra_result_reuse(void)
{
    const double adjacency[] = {
        0.0, 2.0, 5.0,
        2.0, 0.0, 1.0,
        5.0, 1.0, 0.0
    };
    WeightedGraph graph = {3, (double *)adjacency};
    ShortestPathResult result = {0, NULL, NULL};
    double *original_distances;
    ptrdiff_t *original_predecessors;
    size_t i;
    int valid;

    valid = greedy_dijkstra_matrix(&graph, 0, &result);
    if (!valid) return 0;
    original_distances = result.distances;
    original_predecessors = result.predecessors;
    valid = !greedy_dijkstra_matrix(&graph, 3, &result) &&
            result.vertex_count == 3 &&
            result.distances == original_distances &&
            result.predecessors == original_predecessors &&
            result.distances[2] == 3.0 &&
            result.predecessors[2] == 1;
    for (i = 0; valid && i < 32; ++i) {
        valid = greedy_dijkstra_matrix(&graph, 0, &result) &&
                result.vertex_count == 3 &&
                result.distances[2] == 3.0 &&
                result.predecessors[2] == 1;
    }
    greedy_shortest_path_destroy(&result);
    return valid && result.distances == NULL &&
           result.predecessors == NULL && result.vertex_count == 0;
}

static int test_greedy_prim(void)
{
    double adjacency[16] = {
        0.0, 1.0, 5.0, DBL_MAX,
        1.0, 0.0, 2.0, 3.0,
        5.0, 2.0, 0.0, 4.0,
        DBL_MAX, 3.0, 4.0, 0.0
    };
    WeightedGraph graph = {4, adjacency};
    PrimResult result = {0, NULL, 0.0};

    if (!greedy_prim_matrix(&graph, 0, &result) ||
        result.edge_count != 3 ||
        fabs(result.total_weight - 6.0) > 1e-9) {
        greedy_prim_destroy(&result);
        return 0;
    }
    greedy_prim_destroy(&result);
    return 1;
}

static int test_prim_result_reuse(void)
{
    const double adjacency[] = {
        0.0, 2.0, 5.0,
        2.0, 0.0, 1.0,
        5.0, 1.0, 0.0
    };
    WeightedGraph graph = {3, (double *)adjacency};
    PrimResult result = {0, NULL, 0.0};
    WeightedEdge *original_edges;
    size_t i;
    int valid = greedy_prim_matrix(&graph, 0, &result) &&
                result.edge_count == 2 && result.total_weight == 3.0;

    if (!valid) {
        greedy_prim_destroy(&result);
        return 0;
    }
    original_edges = result.edges;
    valid = !greedy_prim_matrix(&graph, 3, &result) &&
            result.edge_count == 2 && result.edges == original_edges &&
            result.total_weight == 3.0 &&
            result.edges[0].weight + result.edges[1].weight == 3.0;
    for (i = 0; valid && i < 32; ++i) {
        valid = greedy_prim_matrix(&graph, 0, &result) &&
                result.edge_count == 2 && result.total_weight == 3.0;
    }
    greedy_prim_destroy(&result);
    return valid && result.edges == NULL && result.edge_count == 0;
}

static int test_greedy_validation(void)
{
    double adjacency[4] = {0.0, -1.0, -1.0, 0.0};
    WeightedGraph graph = {2, adjacency};
    ShortestPathResult path = {0, NULL, NULL};
    PrimResult mst = {0, NULL, 0.0};
    int valid = !greedy_dijkstra_matrix(&graph, 0, &path) &&
                greedy_prim_matrix(&graph, 0, &mst) &&
                mst.edge_count == 1 && mst.total_weight == -1.0;

    greedy_shortest_path_destroy(&path);
    greedy_prim_destroy(&mst);
    graph.vertex_count = SIZE_MAX;
    valid = valid && !greedy_dijkstra_matrix(&graph, 0, &path) &&
            !greedy_prim_matrix(&graph, 0, &mst);
    return valid;
}

static int test_prim_negative_weights(void)
{
    const double adjacency[] = {
        0.0, -5.0, 1.0,
        -5.0, 0.0, -2.0,
        1.0, -2.0, 0.0
    };
    WeightedGraph graph = {3, (double *)adjacency};
    PrimResult result = {0, NULL, 0.0};
    int valid = greedy_prim_matrix(&graph, 0, &result) &&
                result.edge_count == 2 && result.total_weight == -7.0;

    greedy_prim_destroy(&result);
    return valid;
}

static int test_prim_zero_weights(void)
{
    const double adjacency[] = {
        0.0, 0.0, 2.0,
        0.0, 0.0, 0.0,
        2.0, 0.0, 0.0
    };
    WeightedGraph graph = {3, (double *)adjacency};
    PrimResult result = {0, NULL, 0.0};
    int valid = greedy_prim_matrix(&graph, 0, &result) &&
                result.edge_count == 2 && result.total_weight == 0.0;

    greedy_prim_destroy(&result);
    return valid;
}

static int test_prim_disconnected_graph(void)
{
    const double adjacency[] = {
        0.0, -2.0, DBL_MAX, DBL_MAX,
        -2.0, 0.0, DBL_MAX, DBL_MAX,
        DBL_MAX, DBL_MAX, 0.0, 0.0,
        DBL_MAX, DBL_MAX, 0.0, 0.0
    };
    WeightedGraph graph = {4, (double *)adjacency};
    PrimResult result = {0, NULL, 0.0};
    int valid = greedy_prim_matrix(&graph, 0, &result) &&
                result.edge_count == 2 && result.total_weight == -2.0;

    greedy_prim_destroy(&result);
    return valid;
}

static int test_prim_invalid_matrix(void)
{
    double asymmetric[] = {0.0, 1.0, 2.0, 0.0};
    double not_a_number[] = {0.0, NAN, NAN, 0.0};
    double infinite[] = {0.0, INFINITY, INFINITY, 0.0};
    WeightedGraph graph = {2, asymmetric};
    PrimResult result = {0, NULL, 0.0};
    int valid = !greedy_prim_matrix(&graph, 0, &result);

    graph.adjacency = not_a_number;
    valid = valid && !greedy_prim_matrix(&graph, 0, &result);
    graph.adjacency = infinite;
    valid = valid && !greedy_prim_matrix(&graph, 0, &result);
    greedy_prim_destroy(&result);
    return valid;
}

static int test_kruskal_connected_graph(void)
{
    const WeightedEdge edges[] = {
        {0, 1, 1.0}, {0, 2, 4.0}, {1, 2, 2.0},
        {1, 3, 5.0}, {2, 3, 1.0}
    };
    KruskalResult result = {0, NULL, 0.0, 0};
    int valid = greedy_kruskal(4, edges, 5, &result) &&
                result.connected && result.edge_count == 3 &&
                fabs(result.total_weight - 4.0) < 1e-9;
    greedy_kruskal_destroy(&result);
    return valid;
}

static int test_kruskal_disconnected_graph(void)
{
    const WeightedEdge edges[] = {{0, 1, 2.0}, {2, 3, 1.0}};
    KruskalResult result = {0, NULL, 0.0, 0};
    int valid = greedy_kruskal(5, edges, 2, &result) &&
                !result.connected && result.edge_count == 2 &&
                fabs(result.total_weight - 3.0) < 1e-9;
    greedy_kruskal_destroy(&result);
    return valid;
}

static int test_kruskal_negative_weights(void)
{
    const WeightedEdge edges[] = {
        {0, 1, -3.0}, {1, 2, -2.0}, {0, 2, 1.0}
    };
    KruskalResult result = {0, NULL, 0.0, 0};
    int valid = greedy_kruskal(3, edges, 3, &result) &&
                result.connected && result.edge_count == 2 &&
                fabs(result.total_weight + 5.0) < 1e-9;
    greedy_kruskal_destroy(&result);
    return valid;
}

static int test_kruskal_duplicate_edges(void)
{
    const WeightedEdge edges[] = {
        {0, 1, 2.0}, {0, 1, 2.0}, {1, 0, 2.0},
        {1, 2, 1.0}, {0, 2, 4.0}
    };
    KruskalResult result = {0, NULL, 0.0, 0};
    int valid = greedy_kruskal(3, edges, 5, &result) &&
                result.connected && result.edge_count == 2 &&
                fabs(result.total_weight - 3.0) < 1e-9;
    greedy_kruskal_destroy(&result);
    return valid;
}

static int test_kruskal_self_loops(void)
{
    const WeightedEdge edges[] = {
        {0, 0, -100.0}, {1, 1, -50.0}, {0, 1, 3.0}
    };
    KruskalResult result = {0, NULL, 0.0, 0};
    int valid = greedy_kruskal(2, edges, 3, &result) &&
                result.connected && result.edge_count == 1 &&
                result.edges[0].from == 0 && result.edges[0].to == 1 &&
                result.total_weight == 3.0;
    greedy_kruskal_destroy(&result);
    return valid;
}

static int test_kruskal_single_vertex(void)
{
    const WeightedEdge loop[] = {{0, 0, -1.0}};
    KruskalResult result = {0, NULL, 0.0, 0};
    int valid = greedy_kruskal(1, loop, 1, &result) &&
                result.connected && result.edge_count == 0 &&
                result.edges == NULL && result.total_weight == 0.0;
    greedy_kruskal_destroy(&result);
    return valid;
}

static int test_kruskal_empty_graph(void)
{
    KruskalResult result = {0, NULL, 0.0, 0};
    int valid = greedy_kruskal(0, NULL, 0, &result) &&
                !result.connected && result.edge_count == 0 &&
                result.edges == NULL && result.total_weight == 0.0;
    greedy_kruskal_destroy(&result);
    return valid;
}

static int test_kruskal_invalid_input(void)
{
    const WeightedEdge bad_endpoint[] = {{0, 2, 1.0}};
    const WeightedEdge non_finite[] = {{0, 1, INFINITY}};
    KruskalResult result = {0, NULL, 0.0, 0};

    return !greedy_kruskal(2, NULL, 1, &result) &&
           !greedy_kruskal(2, bad_endpoint, 1, &result) &&
           !greedy_kruskal(2, non_finite, 1, &result) &&
           !greedy_kruskal(0, bad_endpoint, 1, &result) &&
           !greedy_kruskal(SIZE_MAX, NULL, 0, &result) &&
           !greedy_kruskal(2, NULL, 0, NULL);
}

static int test_kruskal_result_reuse(void)
{
    const WeightedEdge edges[] = {
        {0, 1, 2.0}, {1, 2, 1.0}, {0, 2, 5.0}
    };
    const WeightedEdge bad_endpoint[] = {{0, 3, 1.0}};
    KruskalResult result = {0, NULL, 0.0, 0};
    WeightedEdge *original_edges;
    size_t i;
    int valid = greedy_kruskal(3, edges, 3, &result) &&
                result.connected && result.edge_count == 2 &&
                result.total_weight == 3.0;

    if (!valid) {
        greedy_kruskal_destroy(&result);
        return 0;
    }
    original_edges = result.edges;
    valid = !greedy_kruskal(3, bad_endpoint, 1, &result) &&
            result.connected && result.edge_count == 2 &&
            result.edges == original_edges && result.total_weight == 3.0 &&
            result.edges[0].weight + result.edges[1].weight == 3.0;
    for (i = 0; valid && i < 32; ++i) {
        valid = greedy_kruskal(3, edges, 3, &result) &&
                result.connected && result.edge_count == 2 &&
                result.total_weight == 3.0;
    }
    greedy_kruskal_destroy(&result);
    return valid && result.edges == NULL && result.edge_count == 0 &&
           !result.connected;
}

static int compare_prim_kruskal_weights(size_t vertex_count,
                                       const WeightedEdge *edges,
                                       size_t edge_count,
                                       const double *adjacency)
{
    WeightedGraph graph;
    PrimResult prim = {0, NULL, 0.0};
    KruskalResult kruskal = {0, NULL, 0.0, 0};
    int valid;

    graph.vertex_count = vertex_count;
    graph.adjacency = (double *)adjacency;
    valid = greedy_prim_matrix(&graph, 0, &prim) &&
            greedy_kruskal(vertex_count, edges, edge_count, &kruskal) &&
            kruskal.connected &&
            fabs(prim.total_weight - kruskal.total_weight) < 1e-9;
    greedy_prim_destroy(&prim);
    greedy_kruskal_destroy(&kruskal);
    return valid;
}

static int test_kruskal_matches_prim(void)
{
    const WeightedEdge first_edges[] = {
        {0, 1, 1.0}, {0, 2, 4.0}, {1, 2, 2.0},
        {1, 3, 5.0}, {2, 3, 1.0}
    };
    const double first_matrix[] = {
        0.0, 1.0, 4.0, DBL_MAX,
        1.0, 0.0, 2.0, 5.0,
        4.0, 2.0, 0.0, 1.0,
        DBL_MAX, 5.0, 1.0, 0.0
    };
    const WeightedEdge second_edges[] = {
        {0, 1, 7.0}, {0, 2, 3.0}, {0, 3, 8.0}, {0, 4, 5.0},
        {1, 2, 4.0}, {1, 3, 2.0}, {1, 4, 6.0},
        {2, 3, 9.0}, {2, 4, 1.0}, {3, 4, 3.0}
    };
    const double second_matrix[] = {
        0.0, 7.0, 3.0, 8.0, 5.0,
        7.0, 0.0, 4.0, 2.0, 6.0,
        3.0, 4.0, 0.0, 9.0, 1.0,
        8.0, 2.0, 9.0, 0.0, 3.0,
        5.0, 6.0, 1.0, 3.0, 0.0
    };

    return compare_prim_kruskal_weights(4, first_edges, 5, first_matrix) &&
           compare_prim_kruskal_weights(5, second_edges, 10, second_matrix);
}

static int test_prim_kruskal_negative_weights(void)
{
    const WeightedEdge first_edges[] = {
        {0, 1, -4.0}, {0, 2, 2.0}, {1, 2, -1.0}
    };
    const double first_matrix[] = {
        0.0, -4.0, 2.0,
        -4.0, 0.0, -1.0,
        2.0, -1.0, 0.0
    };
    const WeightedEdge second_edges[] = {
        {0, 1, -2.0}, {0, 2, 3.0}, {0, 3, -1.0},
        {1, 2, -3.0}, {1, 3, 4.0}, {2, 3, -2.0}
    };
    const double second_matrix[] = {
        0.0, -2.0, 3.0, -1.0,
        -2.0, 0.0, -3.0, 4.0,
        3.0, -3.0, 0.0, -2.0,
        -1.0, 4.0, -2.0, 0.0
    };

    return compare_prim_kruskal_weights(3, first_edges, 3, first_matrix) &&
           compare_prim_kruskal_weights(4, second_edges, 6, second_matrix);
}

static int string_match_result_equals(const StringMatchResult *result,
                                      const size_t *positions,
                                      size_t count)
{
    return result->count == count &&
           result->capacity >= result->count &&
           ((count == 0) == (result->positions == NULL)) &&
           (count == 0 ||
            memcmp(result->positions, positions,
                   count * sizeof(*positions)) == 0);
}

static int string_match_all_equal(const char *text,
                                  size_t text_length,
                                  const char *pattern,
                                  size_t pattern_length,
                                  const size_t *expected,
                                  size_t expected_count)
{
    StringMatchResult naive;
    StringMatchResult rabin_karp;
    StringMatchResult kmp;
    StringMatchResult finite_automaton;
    int valid;

    string_match_result_init(&naive);
    string_match_result_init(&rabin_karp);
    string_match_result_init(&kmp);
    string_match_result_init(&finite_automaton);
    valid = string_match_naive(text, text_length, pattern, pattern_length,
                               &naive) &&
            string_match_rabin_karp(text, text_length, pattern,
                                    pattern_length, &rabin_karp) &&
            string_match_kmp(text, text_length, pattern, pattern_length,
                             &kmp) &&
            string_match_finite_automaton(text, text_length, pattern,
                                          pattern_length,
                                          &finite_automaton) &&
            string_match_result_equals(&naive, expected, expected_count) &&
            string_match_result_equals(&rabin_karp, expected,
                                       expected_count) &&
            string_match_result_equals(&kmp, expected, expected_count) &&
            string_match_result_equals(&finite_automaton, expected,
                                       expected_count);
    string_match_result_destroy(&naive);
    string_match_result_destroy(&rabin_karp);
    string_match_result_destroy(&kmp);
    string_match_result_destroy(&finite_automaton);
    return valid;
}

static int test_naive_string_matching(void)
{
    const size_t overlapping[] = {0, 2, 4};
    const size_t at_end[] = {4};
    const size_t repeated[] = {0, 1, 2};
    const size_t empty_pattern[] = {0, 1, 2, 3};
    const size_t empty_text_pattern[] = {0};

    return string_match_all_equal("ABABABA", 7, "ABA", 3,
                                  overlapping, 3) &&
           string_match_all_equal("ABCDE", 5, "E", 1, at_end, 1) &&
           string_match_all_equal("AAAA", 4, "AA", 2, repeated, 3) &&
           string_match_all_equal("ABC", 3, "LONG", 4, NULL, 0) &&
           string_match_all_equal("ABC", 3, "Z", 1, NULL, 0) &&
           string_match_all_equal("ABC", 3, "", 0, empty_pattern, 4) &&
           string_match_all_equal("", 0, "", 0,
                                  empty_text_pattern, 1) &&
           string_match_all_equal("", 0, "A", 1, NULL, 0);
}

static int test_rabin_karp_string_matching(void)
{
    const unsigned char collision_text[] = {0x3B, 0x9A, 0xCA, 0x07};
    const unsigned char zero_pattern[] = {0, 0, 0, 0};
    const size_t rolling_matches[] = {2, 4};
    uint64_t text_hash = 0;
    uint64_t pattern_hash = 0;
    size_t i;

    if (!string_match_all_equal("xxABABA", 7, "ABA", 3,
                                rolling_matches, 2)) {
        return 0;
    }
    for (i = 0; i < sizeof(collision_text); ++i) {
        text_hash = (text_hash * UINT64_C(256) + collision_text[i]) %
                    UINT64_C(1000000007);
        pattern_hash = (pattern_hash * UINT64_C(256) + zero_pattern[i]) %
                       UINT64_C(1000000007);
    }
    if (text_hash != pattern_hash) {
        return 0;
    }

    {
        StringMatchResult result;
        int valid;

        string_match_result_init(&result);
        valid = string_match_rabin_karp(
                    (const char *)collision_text, sizeof(collision_text),
                    (const char *)zero_pattern, sizeof(zero_pattern),
                    &result) &&
                result.count == 0 && result.positions == NULL;
        string_match_result_destroy(&result);
        return valid;
    }
}

static int test_kmp_string_matching(void)
{
    const size_t repeated_prefix[] = {0, 1};
    const size_t repeated_pattern[] = {2, 8};
    const size_t full_pattern[] = {0};

    return string_match_all_equal("AAAAA", 5, "AAAA", 4,
                                  repeated_prefix, 2) &&
           string_match_all_equal("ZZABABACABABAC", 14, "ABABAC", 6,
                                  repeated_pattern, 2) &&
           string_match_all_equal("ABABCABAB", 9, "ABABCABAB", 9,
                                  full_pattern, 1);
}

static int test_finite_automaton_string_matching(void)
{
    const unsigned char text[] = {0, 0xFF, 0, 0xFF, 0};
    const unsigned char pattern[] = {0xFF, 0};
    const size_t binary_matches[] = {1, 3};
    const size_t overlaps[] = {0, 2, 4};

    return string_match_all_equal((const char *)text, sizeof(text),
                                  (const char *)pattern, sizeof(pattern),
                                  binary_matches, 2) &&
           string_match_all_equal("ABABABABA", 9, "ABABA", 5,
                                  overlaps, 3);
}

static int test_string_matching_cross_validation(void)
{
    char text[96];
    char pattern[24];
    char *large_text;
    const char large_pattern[] = "ABABAC";
    const size_t large_match[] = {131066};
    uint32_t state = UINT32_C(0x5EED1234);
    size_t trial;
    StringMatchResult saved;
    size_t *saved_positions;
    int valid = 1;

    for (trial = 0; trial < 128; ++trial) {
        size_t text_length;
        size_t pattern_length;
        size_t expected[97];
        size_t expected_count = 0;
        size_t i;

        state = state * UINT32_C(1664525) + UINT32_C(1013904223);
        text_length = (size_t)(state % 96);
        state = state * UINT32_C(1664525) + UINT32_C(1013904223);
        pattern_length = (size_t)(state % 24);
        for (i = 0; i < text_length; ++i) {
            state = state * UINT32_C(1664525) + UINT32_C(1013904223);
            text[i] = (char)((state >> 24) & UINT32_C(7));
        }
        for (i = 0; i < pattern_length; ++i) {
            state = state * UINT32_C(1664525) + UINT32_C(1013904223);
            pattern[i] = (char)((state >> 24) & UINT32_C(7));
        }

        if (pattern_length == 0) {
            for (i = 0; i <= text_length; ++i) {
                expected[expected_count++] = i;
            }
        } else if (pattern_length <= text_length) {
            size_t start;
            size_t windows = text_length - pattern_length + 1;

            for (start = 0; start < windows; ++start) {
                size_t offset = 0;
                while (offset < pattern_length &&
                       text[start + offset] == pattern[offset]) {
                    ++offset;
                }
                if (offset == pattern_length) {
                    expected[expected_count++] = start;
                }
            }
        }

        if (!string_match_all_equal(text, text_length, pattern,
                                    pattern_length, expected,
                                    expected_count)) {
            return 0;
        }
    }

    large_text = (char *)malloc(131072);
    if (large_text == NULL) {
        return 0;
    }
    memset(large_text, 'A', 131072);
    memcpy(large_text + 131066, large_pattern, sizeof(large_pattern) - 1);
    valid = string_match_all_equal(large_text, 131072, large_pattern,
                                   sizeof(large_pattern) - 1,
                                   large_match, 1);
    free(large_text);
    if (!valid) {
        return 0;
    }

    string_match_result_init(&saved);
    valid = string_match_naive("ABC", 3, "B", 1, &saved);
    saved_positions = saved.positions;
    valid = valid &&
            !string_match_naive(NULL, 3, "B", 1, &saved) &&
            !string_match_rabin_karp("ABC", 3, NULL, 1, &saved) &&
            !string_match_kmp("ABC", 3, "B", 1, NULL) &&
            !string_match_finite_automaton("ABC", 3, NULL, 1, &saved) &&
            !string_match_naive("x", SIZE_MAX, "", 0, &saved) &&
            !string_match_kmp("x", SIZE_MAX, "x", SIZE_MAX, &saved) &&
            !string_match_finite_automaton(
                "x", SIZE_MAX, "x", SIZE_MAX, &saved) &&
            saved.positions == saved_positions && saved.count == 1 &&
            saved.positions[0] == 1;
    valid = valid &&
            string_match_kmp("ABC", 3, "Z", 1, &saved) &&
            saved.count == 0 && saved.positions == NULL;
    string_match_result_destroy(&saved);
    return valid;
}

static int test_sequence_statistics(void)
{
    const char bases[] = "ACGTN";
    SequenceStatistics statistics;

    if (!sequence_statistics_compute(bases, 5, &statistics) ||
        statistics.a_count != 1 || statistics.c_count != 1 ||
        statistics.g_count != 1 || statistics.t_count != 1 ||
        statistics.ambiguous_count != 1 ||
        statistics.gc_percentage_defined) {
        return 0;
    }
    if (!sequence_statistics_compute("GCGC", 4, &statistics) ||
        !statistics.gc_percentage_defined ||
        statistics.gc_percentage != 100.0) {
        return 0;
    }
    return sequence_statistics_compute(NULL, 0, &statistics) &&
           statistics.a_count == 0 && !statistics.gc_percentage_defined &&
           !sequence_statistics_compute(NULL, 1, &statistics);
}

typedef struct {
    size_t current_runs;
    size_t total_runs;
} BenchmarkTestContext;

static int benchmark_test_prepare(void *context_pointer)
{
    BenchmarkTestContext *context =
        (BenchmarkTestContext *)context_pointer;
    context->current_runs = 0;
    return 1;
}

static int benchmark_test_operation(void *context_pointer)
{
    BenchmarkTestContext *context =
        (BenchmarkTestContext *)context_pointer;
    ++context->current_runs;
    ++context->total_runs;
    return 1;
}

static int benchmark_test_verify(void *context_pointer)
{
    BenchmarkTestContext *context =
        (BenchmarkTestContext *)context_pointer;
    return context->current_runs == 1;
}

static void benchmark_test_cleanup(void *context_pointer)
{
    (void)context_pointer;
}

static int test_benchmark_measurement(void)
{
    BenchmarkTestContext context = {0, 0};
    BenchmarkRecord record;

    if (!benchmark_measure("unit callback", 12, 3,
                           benchmark_test_prepare,
                           benchmark_test_operation,
                           benchmark_test_verify,
                           benchmark_test_cleanup, NULL, &context, &record) ||
        context.total_runs != 3 || record.repetitions != 3 ||
        record.input_size != 12 || !record.verified ||
        record.average_seconds < 0.0 ||
        benchmark_measure("", 1, 1, benchmark_test_prepare,
                          benchmark_test_operation, benchmark_test_verify,
                          benchmark_test_cleanup, NULL, &context, &record)) {
        return 0;
    }
    return 1;
}

static int test_benchmark_suite(void)
{
    static const char *datasets[] = {
        "small", "medium", "large", "adversarial"
    };
    static const char *algorithms[] = {
        "Naive", "Rabin-Karp", "KMP", "Finite Automaton"
    };
    BenchmarkRecord records[BENCHMARK_SUITE_MAX_RECORDS];
    size_t count = 0;
    size_t i;

    if (!benchmark_run_suite(records, BENCHMARK_SUITE_MAX_RECORDS, &count) ||
        count != BENCHMARK_SUITE_MAX_RECORDS) {
        return 0;
    }
    for (i = 0; i < count; ++i) {
        if (records[i].algorithm[0] == '\0' ||
            records[i].average_seconds < 0.0 ||
            records[i].repetitions == 0 || !records[i].verified) {
            return 0;
        }
    }
    for (i = 0; i < sizeof(datasets) / sizeof(datasets[0]); ++i) {
        size_t algorithm_index;
        size_t first_record = 8 + i * 4;

        for (algorithm_index = 0; algorithm_index < 4; ++algorithm_index) {
            BenchmarkRecord *record = &records[first_record + algorithm_index];

            if (strcmp(record->dataset, datasets[i]) != 0 ||
                strcmp(record->algorithm, algorithms[algorithm_index]) != 0 ||
                record->input_size != records[first_record].input_size ||
                !record->has_measured_result_bytes ||
                !record->has_estimated_auxiliary_bytes ||
                !record->has_result_count) {
                return 0;
            }
        }
    }
    return !benchmark_run_suite(records, BENCHMARK_SUITE_MAX_RECORDS - 1,
                                &count);
}

static int file_contains_text(const char *path, const char *needle)
{
    FILE *file = fopen(path, "rb");
    char content[4096];
    size_t bytes_read;

    if (file == NULL) {
        return 0;
    }
    bytes_read = fread(content, 1, sizeof(content) - 1, file);
    content[bytes_read] = '\0';
    fclose(file);
    return strstr(content, needle) != NULL;
}

static int test_result_exports(void)
{
    char identifier[] = "id,\"quoted\"";
    char bases[] = "ACGNT";
    DnaSequence sequence = {identifier, bases, 5};
    BenchmarkRecord benchmark = {
        .algorithm = "test, \"run\"",
        .dataset = "small",
        .input_size = 4,
        .average_seconds = 0.25,
        .repetitions = 2,
        .verified = 1,
        .operation_count = 7,
        .has_operation_count = 1,
        .measured_result_bytes = 16,
        .has_measured_result_bytes = 1,
        .estimated_auxiliary_bytes = 32,
        .has_estimated_auxiliary_bytes = 1,
        .result_count = 2,
        .has_result_count = 1,
        .solution_value = -3,
        .has_solution_value = 1
    };
    AlgorithmComparisonRecord comparison = {"KMP, \"verified\"", 2, 1};
    int successful = 0;

    if (!export_benchmark_records_csv("test-m8-benchmark.csv", &benchmark,
                                      1) ||
        !export_benchmark_records_json("test-m8-benchmark.json", &benchmark,
                                       1) ||
        !export_sequence_statistics_csv("test-m8-sequences.csv", &sequence,
                                        1) ||
        !export_sequence_statistics_json("test-m8-sequences.json", &sequence,
                                         1) ||
        !export_algorithm_comparisons_csv("test-m8-comparisons.csv",
                                          &comparison, 1) ||
        !export_algorithm_comparisons_json("test-m8-comparisons.json",
                                           &comparison, 1)) {
        goto cleanup;
    }
    successful =
        file_contains_text("test-m8-benchmark.csv", "\"test, \"\"run\"\"\"") &&
        file_contains_text("test-m8-benchmark.csv",
                           "measured_result_bytes,estimated_auxiliary_bytes") &&
        file_contains_text("test-m8-benchmark.csv", ",16,32,2,-3") &&
        file_contains_text("test-m8-benchmark.json", "test, \\\"run\\\"") &&
        file_contains_text("test-m8-benchmark.json",
                           "\"measured_result_bytes\":16") &&
        file_contains_text("test-m8-benchmark.json",
                           "\"solution_value\":-3") &&
        file_contains_text("test-m8-sequences.csv", "\"id,\"\"quoted\"\"\"") &&
        file_contains_text("test-m8-sequences.json", "id,\\\"quoted\\\"") &&
        file_contains_text("test-m8-comparisons.csv",
                           "\"KMP, \"\"verified\"\"\"") &&
        file_contains_text("test-m8-comparisons.json",
                           "KMP, \\\"verified\\\"");

cleanup:
    (void)remove("test-m8-benchmark.csv");
    (void)remove("test-m8-benchmark.json");
    (void)remove("test-m8-sequences.csv");
    (void)remove("test-m8-sequences.json");
    (void)remove("test-m8-comparisons.csv");
    (void)remove("test-m8-comparisons.json");
    return successful;
}

static int test_cli_menu_integration(void)
{
    static const char input_text[] =
        "99\n"
        "01\n99\n00\n"
        "00\nn\n"
        "01\n01\nACGT\n00\n"
        "02\n01\n\n00\n"
        "03\n01\n\n1\n9\n1\n00\n"
        "04\n01\n00\n"
        "05\n01\n00\n"
        "06\n01\nACGT\nAGT\n00\n"
        "07\n01\n15\n00\n"
        "08\n01\n00\n"
        "09\n01\n\n1\nACG\n4\n00\n"
        "10\n01\n00\n"
        "11\n01\n\nmenu-test-m8\n00\n"
        "12\n03\n04\n01\nACG\nAG\n\n\n\n"
        "02\nA\nA\n3\n-2\n-4\n04\n00\n"
        "00\ny\n";
    FILE *input = tmpfile();
    FILE *output = tmpfile();
    char transcript[65536];
    size_t output_length;
    int successful;

    if (input == NULL || output == NULL) {
        if (input != NULL) {
            fclose(input);
        }
        if (output != NULL) {
            fclose(output);
        }
        return 0;
    }
    fputs(input_text, input);
    rewind(input);
    successful = genome_cli_run(input, output);
    rewind(output);
    output_length = fread(transcript, 1, sizeof(transcript) - 1, output);
    transcript[output_length] = '\0';
    successful = successful &&
                 strstr(transcript, "PROJECT MODULES") != NULL &&
                 strstr(transcript,
                        "[10] Benchmarking and Performance Analysis") != NULL &&
                 strstr(transcript, "[11] Export Results (CSV / JSON)") !=
                     NULL &&
                 strstr(transcript, "[12] DNA Sequence Alignment") != NULL &&
                 strstr(transcript, "Invalid menu choice") != NULL &&
                 strstr(transcript, "Exit cancelled.") != NULL &&
                 strstr(transcript, "A=1 C=1 G=1 T=1") != NULL &&
                 strstr(transcript, "Found sample_beta at sorted index") != NULL &&
                 strstr(transcript, "LCS length=3 sequence=AGT") != NULL &&
                 strstr(transcript, "Algorithm       : Needleman-Wunsch") !=
                     NULL &&
                 strstr(transcript, "Alignment Score : 2") != NULL &&
                 strstr(transcript, "Gaps            : 1") != NULL &&
                 strstr(transcript, "Algorithm       : Smith-Waterman") !=
                     NULL &&
                 strstr(transcript, "Alignment Score : 3") != NULL &&
                 strstr(transcript, "ALIGNMENT HISTORY (THIS SESSION)") !=
                     NULL &&
                 strstr(transcript, "Needleman-Wunsch score=2") != NULL &&
                 strstr(transcript,
                        "Enter match reward (points for identical DNA "
                        "characters) [+2]") != NULL &&
                 strstr(transcript, "Textbook flow network maximum flow: 23") !=
                     NULL &&
                 strstr(transcript,
                        "Maximum bipartite matching size: 3.") != NULL &&
                 strstr(transcript, "All algorithms agree: yes.") != NULL &&
                 strstr(transcript, "Compare exact 0/1 Knapsack search") !=
                     NULL &&
                 strstr(transcript, "menu-test-m8-comparisons.json") != NULL &&
                 strstr(transcript, "Goodbye.") != NULL;
    (void)remove("menu-test-m8-sequences.csv");
    (void)remove("menu-test-m8-sequences.json");
    (void)remove("menu-test-m8-benchmarks.csv");
    (void)remove("menu-test-m8-benchmarks.json");
    (void)remove("menu-test-m8-comparisons.csv");
    (void)remove("menu-test-m8-comparisons.json");
    fclose(input);
    fclose(output);
    return successful;
}

static int demo_output_contains(const char *scenario,
                                const char *expected,
                                const char *cleanup_prefix)
{
    FILE *output = tmpfile();
    char transcript[8192];
    size_t length;
    int successful;

    if (output == NULL) {
        return 0;
    }
    successful = genome_cli_run_demo(output, scenario);
    rewind(output);
    length = fread(transcript, 1, sizeof(transcript) - 1, output);
    transcript[length] = '\0';
    successful = successful && strstr(transcript, expected) != NULL;
    fclose(output);
    if (cleanup_prefix != NULL) {
        char path[128];
        const char *suffixes[] = {
            "sequences.csv", "sequences.json", "matches.csv", "matches.json",
            "benchmarks.csv", "benchmarks.json"
        };
        size_t i;
        for (i = 0; i < sizeof(suffixes) / sizeof(suffixes[0]); ++i) {
            (void)snprintf(path, sizeof(path), "%s-%s", cleanup_prefix,
                           suffixes[i]);
            (void)remove(path);
        }
    }
    return successful;
}

static int test_cli_normal_demo(void)
{
    int successful =
        demo_output_contains("normal",
                             "Exported sequence statistics and "
                             "string-match comparison results.",
                             NULL) &&
        file_contains_text("demo-normal-matches.csv", "Naive,2,true") &&
        file_contains_text("demo-normal-matches.csv", "Rabin-Karp,2,true") &&
        file_contains_text("demo-normal-matches.csv", "KMP,2,true") &&
        file_contains_text("demo-normal-matches.csv",
                           "Finite Automaton,2,true") &&
        file_contains_text("demo-normal-matches.json",
                           "\"algorithm\":\"Naive\",\"result_count\":2,"
                           "\"matches_reference\":true") &&
        file_contains_text("demo-normal-matches.json",
                           "\"algorithm\":\"Rabin-Karp\",\"result_count\":2,"
                           "\"matches_reference\":true") &&
        file_contains_text("demo-normal-matches.json",
                           "\"algorithm\":\"KMP\",\"result_count\":2,"
                           "\"matches_reference\":true") &&
        file_contains_text("demo-normal-matches.json",
                           "\"algorithm\":\"Finite Automaton\","
                           "\"result_count\":2,\"matches_reference\":true");

    (void)remove("demo-normal-sequences.csv");
    (void)remove("demo-normal-sequences.json");
    (void)remove("demo-normal-matches.csv");
    (void)remove("demo-normal-matches.json");
    return successful;
}

static int test_cli_constraint_shock_demo(void)
{
    return demo_output_contains("shock", "Capacity 12 -> maximum value",
                                NULL);
}

static int test_cli_stress_demo(void)
{
    return demo_output_contains("stress",
                                "Exported the twelve verified string-search "
                                "benchmark records.",
                                "demo-stress");
}

static int tests_run = 0;

static int run_test(const char *name, TestFunction test)
{
    ++tests_run;
    if (test()) {
        printf("PASS: %s\n", name);
        return 0;
    }

    fprintf(stderr, "FAIL: %s\n", name);
    return 1;
}

int main(void)
{
    int failures = 0;

    failures += run_test("divide-and-conquer maximum subarray",
                         test_maximum_subarray);
    failures += run_test("dynamic programming 0/1 knapsack",
                         test_dp_knapsack_01);
    failures += run_test("0/1 knapsack selection reconstruction",
                         test_dp_knapsack_01_reconstruction);
    failures += run_test("dynamic programming longest common subsequence",
                         test_dp_longest_common_subsequence);
    failures += run_test("LCS input and allocation validation",
                         test_dp_longest_common_subsequence_validation);
    failures += run_test("Needleman-Wunsch global DNA alignment",
                         test_dna_alignment_global);
    failures += run_test("Smith-Waterman local DNA alignment",
                         test_dna_alignment_local);
    failures += run_test("DNA alignment validation and overflow",
                         test_dna_alignment_validation_and_overflow);
    failures += run_test("dynamic programming matrix chain multiplication",
                         test_dp_matrix_chain_multiplication);
    failures += run_test("dynamic programming optimal binary search tree",
                         test_dp_optimal_bst);
    failures += run_test("optimal binary search tree validation",
                         test_dp_optimal_bst_validation);
    failures += run_test("dynamic programming Bellman-Ford shortest paths",
                         test_dp_bellman_ford_shortest_paths);
    failures += run_test("Bellman-Ford negative cycle detection",
                         test_dp_bellman_ford_negative_cycles);
    failures += run_test("Bellman-Ford validation and result reuse",
                         test_dp_bellman_ford_validation_and_reuse);
    failures += run_test("dynamic programming Floyd-Warshall paths",
                         test_dp_floyd_warshall_paths);
    failures += run_test("Floyd-Warshall negative cycle detection",
                         test_dp_floyd_warshall_negative_cycles);
    failures += run_test("Floyd-Warshall validation and result reuse",
                         test_dp_floyd_warshall_validation_and_reuse);
    failures += run_test("Bellman-Ford and Floyd-Warshall agreement",
                         test_dp_floyd_bellman_ford_agreement);
    failures += run_test("dynamic programming TSP bitmask examples",
                         test_dp_tsp_bitmask_examples);
    failures += run_test("TSP no-tour and maximum-size cases",
                         test_dp_tsp_no_tour_and_limits);
    failures += run_test("TSP validation and result reuse",
                         test_dp_tsp_validation_and_reuse);
    failures += run_test("backtracking Sum of Subsets",
                         test_backtracking_sum_of_subsets);
    failures += run_test("backtracking 0/1 Knapsack",
                         test_backtracking_knapsack);
    failures += run_test("FIFO Branch and Bound 0/1 Knapsack",
                         test_fifo_branch_bound_knapsack);
    failures += run_test("LIFO Branch and Bound 0/1 Knapsack",
                         test_lifo_branch_bound_knapsack);
    failures += run_test("Least-Cost Branch and Bound 0/1 Knapsack",
                         test_least_cost_branch_bound_knapsack);
    failures += run_test("Knapsack algorithm agreement",
                         test_knapsack_algorithm_agreement);
    failures += run_test("Ford-Fulkerson textbook and flow cases",
                         test_ford_fulkerson_flow_examples);
    failures += run_test("Ford-Fulkerson reverse residual adjustment",
                         test_ford_fulkerson_reverse_residual_edges);
    failures += run_test("Ford-Fulkerson validation and overflow",
                         test_ford_fulkerson_validation_and_overflow);
    failures += run_test("maximum bipartite matching examples",
                         test_bipartite_matching_examples);
    failures += run_test("bipartite matching validation",
                         test_bipartite_matching_validation);
    failures += run_test("bipartite matching brute-force agreement",
                         test_bipartite_matching_brute_force_agreement);
    failures += run_test("naive string matching",
                         test_naive_string_matching);
    failures += run_test("Rabin-Karp collision and rolling hash",
                         test_rabin_karp_string_matching);
    failures += run_test("KMP prefix and overlap cases",
                         test_kmp_string_matching);
    failures += run_test("finite-automaton byte matching",
                         test_finite_automaton_string_matching);
    failures += run_test("string matching cross-validation",
                         test_string_matching_cross_validation);
    failures += run_test("exponentiation by squaring",
                         test_fast_exponentiation);
    failures += run_test("celebrity identification", test_celebrity_problem);
    failures += run_test("classical rectangular matrix multiplication",
                         test_classical_matrix_multiplication);
    failures += run_test("Strassen matrix multiplication",
                         test_strassen_matrix_multiplication);
    failures += run_test("merge sort by sequence length",
                         test_merge_sort_by_length);
    failures += run_test("merge sort by identifier",
                         test_merge_sort_by_identifier);
    failures += run_test("quick sort by sequence length",
                         test_quick_sort_by_length);
    failures += run_test("quick sort by identifier",
                         test_quick_sort_by_identifier);
    failures += run_test("iterative binary search by length",
                         test_iterative_binary_search_by_length);
    failures += run_test("iterative binary search by identifier",
                         test_iterative_binary_search_by_identifier);
    failures += run_test("recursive binary search by length",
                         test_recursive_binary_search_by_length);
    failures += run_test("recursive binary search by identifier",
                         test_recursive_binary_search_by_identifier);
    failures += run_test("valid sequences", test_valid_sequences);
    failures += run_test("uppercase sequence", test_uppercase_sequence);
    failures += run_test("lowercase sequence", test_lowercase_sequence);
    failures += run_test("single nucleotide", test_single_nucleotide);
    failures += run_test("invalid DNA characters", test_invalid_characters);
    failures += run_test("empty DNA sequence", test_empty_sequence);
    failures += run_test("NULL DNA input", test_null_sequence);
    failures += run_test("sequence creation and destruction",
                         test_create_and_destroy);
    failures += run_test("creation into occupied structure",
                         test_create_into_occupied_structure);
    failures += run_test("invalid identifier", test_invalid_identifier);
    failures += run_test("invalid sequence creation",
                         test_invalid_sequence_creation);
    failures += run_test("FASTA single record", test_fasta_single_record);
    failures += run_test("FASTA multiple records", test_fasta_multiple_records);
    failures += run_test("FASTA multiline sequence",
                         test_fasta_multiline_sequence);
    failures += run_test("FASTA invalid DNA", test_fasta_invalid_dna);
    failures += run_test("FASTA data before header",
                         test_fasta_missing_header);
    failures += run_test("FASTA empty identifier",
                         test_fasta_empty_identifier);
    failures += run_test("FASTA whitespace-only identifier",
                         test_fasta_whitespace_identifier);
    failures += run_test("FASTA embedded NULL byte",
                         test_fasta_embedded_null);
    failures += run_test("FASTA empty sequence", test_fasta_empty_sequence);
    failures += run_test("FASTA empty file", test_fasta_empty_file);
    failures += run_test("FASTA missing file", test_fasta_missing_file);
    failures += run_test("FASTA NULL output parameters",
                         test_fasta_null_output_parameters);
    failures += run_test("FASTA CRLF line endings", test_fasta_crlf);
    failures += run_test("FASTA long sequence", test_fasta_long_sequence);
    failures += run_test("greedy activity selection",
                         test_greedy_activity_selection);
    failures += run_test("greedy fractional knapsack",
                         test_greedy_fractional_knapsack);
    failures += run_test("greedy Huffman coding", test_greedy_huffman);
    failures += run_test("greedy Huffman input validation",
                         test_greedy_huffman_validation);
    failures += run_test("greedy Huffman edge cases",
                         test_greedy_huffman_edge_cases);
    failures += run_test("greedy Dijkstra shortest paths",
                         test_greedy_dijkstra);
    failures += run_test("Dijkstra result reuse and failure preservation",
                         test_dijkstra_result_reuse);
    failures += run_test("greedy Prim minimum spanning tree",
                         test_greedy_prim);
    failures += run_test("Prim result reuse and failure preservation",
                         test_prim_result_reuse);
    failures += run_test("Prim negative edge weights",
                         test_prim_negative_weights);
    failures += run_test("Prim zero-weight edges", test_prim_zero_weights);
    failures += run_test("Prim disconnected graph",
                         test_prim_disconnected_graph);
    failures += run_test("Prim invalid adjacency matrices",
                         test_prim_invalid_matrix);
    failures += run_test("greedy invalid input validation",
                         test_greedy_validation);
    failures += run_test("Kruskal connected graph",
                         test_kruskal_connected_graph);
    failures += run_test("Kruskal disconnected graph",
                         test_kruskal_disconnected_graph);
    failures += run_test("Kruskal negative edge weights",
                         test_kruskal_negative_weights);
    failures += run_test("Kruskal duplicate edges",
                         test_kruskal_duplicate_edges);
    failures += run_test("Kruskal self-loops", test_kruskal_self_loops);
    failures += run_test("Kruskal single vertex", test_kruskal_single_vertex);
    failures += run_test("Kruskal empty graph", test_kruskal_empty_graph);
    failures += run_test("Kruskal invalid input", test_kruskal_invalid_input);
    failures += run_test("Kruskal result reuse and failure preservation",
                         test_kruskal_result_reuse);
    failures += run_test("Kruskal and Prim weight agreement",
                         test_kruskal_matches_prim);
    failures += run_test("Prim and Kruskal negative-weight agreement",
                         test_prim_kruskal_negative_weights);
    failures += run_test("sequence statistics and GC percentage",
                         test_sequence_statistics);
    failures += run_test("benchmark measurement and verification",
                         test_benchmark_measurement);
    failures += run_test("verified multi-algorithm benchmark suite",
                         test_benchmark_suite);
    failures += run_test("CSV and JSON result exports",
                         test_result_exports);
    failures += run_test("interactive CLI end-to-end workflows",
                         test_cli_menu_integration);
    failures += run_test("normal-operation CLI demonstration",
                         test_cli_normal_demo);
    failures += run_test("constraint-shock CLI demonstration",
                         test_cli_constraint_shock_demo);
    failures += run_test("large-data stress CLI demonstration",
                         test_cli_stress_demo);

    if (failures != 0) {
        fprintf(stderr, "%d of %d test(s) failed.\n", failures, tests_run);
        return 1;
    }

    printf("All %d tests passed.\n", tests_run);
    return 0;
}
