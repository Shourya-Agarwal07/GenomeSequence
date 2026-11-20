#include "genome/algorithms/divide_conquer.h"

static void reset_operations(DcaOperationCounts *operations)
{
    if (operations != NULL) {
        operations->comparisons = 0;
        operations->additions = 0;
        operations->multiplications = 0;
    }
}

static int add_int64_checked(int64_t left, int64_t right, int64_t *sum)
{
    if ((right > 0 && left > INT64_MAX - right) ||
        (right < 0 && left < INT64_MIN - right)) {
        return 0;
    }

    *sum = left + right;
    return 1;
}

static int maximum_subarray_range(const int64_t *values,
                                  size_t start,
                                  size_t end,
                                  MaximumSubarrayResult *result,
                                  DcaOperationCounts *operations)
{
    size_t middle;
    size_t index;
    int64_t running_sum = 0;
    int64_t best_left_sum = INT64_MIN;
    int64_t best_right_sum = INT64_MIN;
    size_t best_left_start = 0;
    size_t best_right_end = 0;
    int64_t crossing_sum;
    MaximumSubarrayResult left_result;
    MaximumSubarrayResult right_result;

    if (end - start == 1) {
        result->sum = values[start];
        result->start_index = start;
        result->end_index = start;
        return 1;
    }

    middle = start + (end - start) / 2;
    if (!maximum_subarray_range(values, start, middle, &left_result,
                                operations) ||
        !maximum_subarray_range(values, middle, end, &right_result,
                                operations)) {
        return 0;
    }

    index = middle;
    while (index > start) {
        --index;
        if (operations != NULL) {
            ++operations->additions;
        }
        if (!add_int64_checked(running_sum, values[index], &running_sum)) {
            return 0;
        }
        if (operations != NULL) {
            ++operations->comparisons;
        }
        if (running_sum > best_left_sum) {
            best_left_sum = running_sum;
            best_left_start = index;
        }
    }

    running_sum = 0;
    for (index = middle; index < end; ++index) {
        if (operations != NULL) {
            ++operations->additions;
        }
        if (!add_int64_checked(running_sum, values[index], &running_sum)) {
            return 0;
        }
        if (operations != NULL) {
            ++operations->comparisons;
        }
        if (running_sum > best_right_sum) {
            best_right_sum = running_sum;
            best_right_end = index;
        }
    }

    if (operations != NULL) {
        ++operations->additions;
    }
    if (!add_int64_checked(best_left_sum, best_right_sum, &crossing_sum)) {
        return 0;
    }

    *result = left_result;
    if (operations != NULL) {
        ++operations->comparisons;
    }
    if (right_result.sum > result->sum) {
        *result = right_result;
    }
    if (operations != NULL) {
        ++operations->comparisons;
    }
    if (crossing_sum > result->sum) {
        result->sum = crossing_sum;
        result->start_index = best_left_start;
        result->end_index = best_right_end;
    }

    return 1;
}

int dca_maximum_subarray(const int64_t *values,
                         size_t count,
                         MaximumSubarrayResult *result,
                         DcaOperationCounts *operations)
{
    reset_operations(operations);
    if (result != NULL) {
        result->sum = 0;
        result->start_index = 0;
        result->end_index = 0;
    }
    if (values == NULL || count == 0 || result == NULL) {
        return 0;
    }

    return maximum_subarray_range(values, 0, count, result, operations);
}
