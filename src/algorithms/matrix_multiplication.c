#include "genome/algorithms/divide_conquer.h"

#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#define STRASSEN_CLASSICAL_THRESHOLD 32

static void reset_operations(DcaOperationCounts *operations)
{
    if (operations != NULL) {
        operations->comparisons = 0;
        operations->additions = 0;
        operations->multiplications = 0;
    }
}

static int matrix_element_count(size_t rows,
                                size_t columns,
                                size_t *element_count)
{
    if (rows != 0 && columns > SIZE_MAX / rows) {
        return 0;
    }
    *element_count = rows * columns;
    return 1;
}

static int matrix_is_valid(const DcaMatrix *matrix)
{
    size_t element_count;

    return matrix != NULL &&
           matrix_element_count(matrix->rows, matrix->columns,
                                &element_count) &&
           (element_count == 0 || matrix->values != NULL);
}

static int matrix_values_are_finite(const DcaMatrix *matrix)
{
    size_t element_count;
    size_t i;

    if (!matrix_element_count(matrix->rows, matrix->columns, &element_count)) {
        return 0;
    }
    for (i = 0; i < element_count; ++i) {
        if (!isfinite(matrix->values[i])) {
            return 0;
        }
    }
    return 1;
}

static int matrix_is_empty(const DcaMatrix *matrix)
{
    return matrix != NULL && matrix->rows == 0 && matrix->columns == 0 &&
           matrix->values == NULL;
}

int dca_matrix_create(DcaMatrix *matrix, size_t rows, size_t columns)
{
    size_t element_count;
    double *values = NULL;

    if (!matrix_is_empty(matrix) ||
        !matrix_element_count(rows, columns, &element_count) ||
        element_count > SIZE_MAX / sizeof(*values)) {
        return 0;
    }
    if (element_count != 0) {
        values = (double *)calloc(element_count, sizeof(*values));
        if (values == NULL) {
            return 0;
        }
    }

    matrix->rows = rows;
    matrix->columns = columns;
    matrix->values = values;
    return 1;
}

void dca_matrix_destroy(DcaMatrix *matrix)
{
    if (matrix != NULL) {
        free(matrix->values);
        matrix->rows = 0;
        matrix->columns = 0;
        matrix->values = NULL;
    }
}

static int multiply_dense(const double *left,
                          size_t rows,
                          size_t inner,
                          const double *right,
                          size_t columns,
                          double *result,
                          DcaOperationCounts *operations)
{
    size_t row;
    size_t column;

    for (row = 0; row < rows; ++row) {
        for (column = 0; column < columns; ++column) {
            double sum = 0.0;
            size_t k;

            if (inner != 0) {
                sum = left[row * inner] * right[column];
                if (operations != NULL) {
                    ++operations->multiplications;
                }
                if (!isfinite(sum)) {
                    return 0;
                }
                for (k = 1; k < inner; ++k) {
                    double product = left[row * inner + k] *
                                     right[k * columns + column];

                    if (operations != NULL) {
                        ++operations->multiplications;
                    }
                    if (!isfinite(product)) {
                        return 0;
                    }
                    sum += product;
                    if (operations != NULL) {
                        ++operations->additions;
                    }
                    if (!isfinite(sum)) {
                        return 0;
                    }
                }
            }
            result[row * columns + column] = sum;
        }
    }
    return 1;
}

int dca_matrix_multiply_classical(const DcaMatrix *left,
                                  const DcaMatrix *right,
                                  DcaMatrix *result,
                                  DcaOperationCounts *operations)
{
    DcaMatrix product = {0, 0, NULL};

    reset_operations(operations);
    if (result == left || result == right || !matrix_is_empty(result) ||
        !matrix_is_valid(left) || !matrix_is_valid(right) ||
        left->columns != right->rows ||
        !matrix_values_are_finite(left) ||
        !matrix_values_are_finite(right)) {
        return 0;
    }
    if (!dca_matrix_create(&product, left->rows, right->columns)) {
        return 0;
    }
    if (!multiply_dense(left->values, left->rows, left->columns,
                        right->values, right->columns, product.values,
                        operations)) {
        dca_matrix_destroy(&product);
        return 0;
    }

    *result = product;
    return 1;
}

static void copy_quadrant(const double *source,
                          size_t source_size,
                          size_t source_row,
                          size_t source_column,
                          double *destination,
                          size_t quadrant_size)
{
    size_t row;

    for (row = 0; row < quadrant_size; ++row) {
        size_t column;

        for (column = 0; column < quadrant_size; ++column) {
            destination[row * quadrant_size + column] =
                source[(source_row + row) * source_size +
                       source_column + column];
        }
    }
}

static int combine_quadrants(const double *first,
                             size_t first_size,
                             size_t first_row,
                             size_t first_column,
                             const double *second,
                             size_t second_size,
                             size_t second_row,
                             size_t second_column,
                             double *destination,
                             size_t quadrant_size,
                             int subtract,
                             DcaOperationCounts *operations)
{
    size_t row;

    for (row = 0; row < quadrant_size; ++row) {
        size_t column;

        for (column = 0; column < quadrant_size; ++column) {
            double left_value =
                first[(first_row + row) * first_size + first_column + column];
            double right_value =
                second[(second_row + row) * second_size +
                       second_column + column];
            double value = subtract ? left_value - right_value
                                    : left_value + right_value;

            if (operations != NULL) {
                ++operations->additions;
            }
            if (!isfinite(value)) {
                return 0;
            }
            destination[row * quadrant_size + column] = value;
        }
    }
    return 1;
}

static int combine_values(double left,
                          double right,
                          int subtract,
                          double *result,
                          DcaOperationCounts *operations)
{
    *result = subtract ? left - right : left + right;
    if (operations != NULL) {
        ++operations->additions;
    }
    return isfinite(*result);
}

static int strassen_recursive(const double *left,
                              const double *right,
                              double *result,
                              size_t size,
                              DcaOperationCounts *operations)
{
    size_t half;
    size_t quadrant_count;
    double *workspace;
    double *sum_left;
    double *sum_right;
    double *product_1;
    double *product_2;
    double *product_3;
    double *product_4;
    double *product_5;
    double *product_6;
    double *product_7;
    size_t row;
    int success = 0;

    if (size <= STRASSEN_CLASSICAL_THRESHOLD) {
        return multiply_dense(left, size, size, right, size, result, operations);
    }

    half = size / 2;
    if (!matrix_element_count(half, half, &quadrant_count) ||
        quadrant_count > SIZE_MAX / (9 * sizeof(*workspace))) {
        return 0;
    }
    workspace = (double *)malloc(9 * quadrant_count * sizeof(*workspace));
    if (workspace == NULL) {
        return 0;
    }

    sum_left = workspace;
    sum_right = sum_left + quadrant_count;
    product_1 = sum_right + quadrant_count;
    product_2 = product_1 + quadrant_count;
    product_3 = product_2 + quadrant_count;
    product_4 = product_3 + quadrant_count;
    product_5 = product_4 + quadrant_count;
    product_6 = product_5 + quadrant_count;
    product_7 = product_6 + quadrant_count;

    if (!combine_quadrants(left, size, 0, 0, left, size, half, half,
                           sum_left, half, 0, operations) ||
        !combine_quadrants(right, size, 0, 0, right, size, half, half,
                           sum_right, half, 0, operations) ||
        !strassen_recursive(sum_left, sum_right, product_1, half, operations)) {
        goto cleanup;
    }

    if (!combine_quadrants(left, size, half, 0, left, size, half, half,
                           sum_left, half, 0, operations)) {
        goto cleanup;
    }
    copy_quadrant(right, size, 0, 0, sum_right, half);
    if (!strassen_recursive(sum_left, sum_right, product_2, half, operations)) {
        goto cleanup;
    }

    copy_quadrant(left, size, 0, 0, sum_left, half);
    if (!combine_quadrants(right, size, 0, half, right, size, half, half,
                           sum_right, half, 1, operations) ||
        !strassen_recursive(sum_left, sum_right, product_3, half, operations)) {
        goto cleanup;
    }

    copy_quadrant(left, size, half, half, sum_left, half);
    if (!combine_quadrants(right, size, half, 0, right, size, 0, 0,
                           sum_right, half, 1, operations) ||
        !strassen_recursive(sum_left, sum_right, product_4, half, operations)) {
        goto cleanup;
    }

    if (!combine_quadrants(left, size, 0, 0, left, size, 0, half,
                           sum_left, half, 0, operations)) {
        goto cleanup;
    }
    copy_quadrant(right, size, half, half, sum_right, half);
    if (!strassen_recursive(sum_left, sum_right, product_5, half, operations)) {
        goto cleanup;
    }

    if (!combine_quadrants(left, size, half, 0, left, size, 0, 0,
                           sum_left, half, 1, operations) ||
        !combine_quadrants(right, size, 0, 0, right, size, 0, half,
                           sum_right, half, 0, operations) ||
        !strassen_recursive(sum_left, sum_right, product_6, half, operations)) {
        goto cleanup;
    }

    if (!combine_quadrants(left, size, 0, half, left, size, half, half,
                           sum_left, half, 1, operations) ||
        !combine_quadrants(right, size, half, 0, right, size, half, half,
                           sum_right, half, 0, operations) ||
        !strassen_recursive(sum_left, sum_right, product_7, half, operations)) {
        goto cleanup;
    }

    for (row = 0; row < half; ++row) {
        size_t column;

        for (column = 0; column < half; ++column) {
            size_t offset = row * half + column;
            double value;

            if (!combine_values(product_1[offset], product_4[offset], 0,
                                &value, operations) ||
                !combine_values(value, product_5[offset], 1,
                                &value, operations) ||
                !combine_values(value, product_7[offset], 0,
                                &result[row * size + column], operations) ||
                !combine_values(product_3[offset], product_5[offset], 0,
                                &result[row * size + half + column],
                                operations) ||
                !combine_values(product_2[offset], product_4[offset], 0,
                                &result[(half + row) * size + column],
                                operations) ||
                !combine_values(product_1[offset], product_2[offset], 1,
                                &value, operations) ||
                !combine_values(value, product_3[offset], 0,
                                &value, operations) ||
                !combine_values(value, product_6[offset], 0,
                                &result[(half + row) * size + half + column],
                                operations)) {
                goto cleanup;
            }
        }
    }
    success = 1;

cleanup:
    free(workspace);
    return success;
}

int dca_matrix_multiply_strassen(const DcaMatrix *left,
                                 const DcaMatrix *right,
                                 DcaMatrix *result,
                                 DcaOperationCounts *operations)
{
    DcaMatrix product = {0, 0, NULL};
    double *padded_left = NULL;
    double *padded_right = NULL;
    double *padded_result = NULL;
    size_t padded_size = 1;
    size_t padded_count;
    size_t row;
    int success = 0;

    reset_operations(operations);
    if (result == left || result == right || !matrix_is_empty(result) ||
        !matrix_is_valid(left) || !matrix_is_valid(right) ||
        left->rows != left->columns ||
        right->rows != right->columns ||
        left->rows != right->rows ||
        !matrix_values_are_finite(left) ||
        !matrix_values_are_finite(right)) {
        return 0;
    }
    if (!dca_matrix_create(&product, left->rows, right->columns)) {
        return 0;
    }
    if (left->rows <= STRASSEN_CLASSICAL_THRESHOLD) {
        if (!multiply_dense(left->values, left->rows, left->columns,
                            right->values, right->columns, product.values,
                            operations)) {
            dca_matrix_destroy(&product);
            return 0;
        }
        *result = product;
        return 1;
    }

    while (padded_size < left->rows) {
        if (padded_size > SIZE_MAX / 2) {
            goto cleanup;
        }
        padded_size *= 2;
    }
    if (!matrix_element_count(padded_size, padded_size, &padded_count) ||
        padded_count > SIZE_MAX / sizeof(*padded_left)) {
        goto cleanup;
    }

    padded_left = (double *)calloc(padded_count, sizeof(*padded_left));
    padded_right = (double *)calloc(padded_count, sizeof(*padded_right));
    padded_result = (double *)calloc(padded_count, sizeof(*padded_result));
    if (padded_left == NULL || padded_right == NULL ||
        padded_result == NULL) {
        goto cleanup;
    }

    for (row = 0; row < left->rows; ++row) {
        size_t column;

        for (column = 0; column < left->columns; ++column) {
            padded_left[row * padded_size + column] =
                left->values[row * left->columns + column];
            padded_right[row * padded_size + column] =
                right->values[row * right->columns + column];
        }
    }

    if (!strassen_recursive(padded_left, padded_right, padded_result,
                            padded_size, operations)) {
        goto cleanup;
    }
    for (row = 0; row < product.rows; ++row) {
        size_t column;

        for (column = 0; column < product.columns; ++column) {
            product.values[row * product.columns + column] =
                padded_result[row * padded_size + column];
        }
    }

    *result = product;
    product.values = NULL;
    success = 1;

cleanup:
    free(padded_left);
    free(padded_right);
    free(padded_result);
    dca_matrix_destroy(&product);
    return success;
}
