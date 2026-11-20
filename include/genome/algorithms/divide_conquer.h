#ifndef GENOME_ALGORITHMS_DIVIDE_CONQUER_H
#define GENOME_ALGORITHMS_DIVIDE_CONQUER_H

#include <stddef.h>
#include <stdint.h>

typedef struct {
    size_t comparisons;
    size_t additions;
    size_t multiplications;
} DcaOperationCounts;

typedef struct {
    int64_t sum;
    size_t start_index;
    size_t end_index;
} MaximumSubarrayResult;

typedef struct {
    size_t rows;
    size_t columns;
    double *values;
} DcaMatrix;

int dca_maximum_subarray(const int64_t *values,
                         size_t count,
                         MaximumSubarrayResult *result,
                         DcaOperationCounts *operations);

int dca_fast_power(double base,
                   int exponent,
                   double *result,
                   DcaOperationCounts *operations);

ptrdiff_t dca_find_celebrity(const unsigned char *relationships,
                            size_t count,
                            DcaOperationCounts *operations);

int dca_matrix_create(DcaMatrix *matrix, size_t rows, size_t columns);
void dca_matrix_destroy(DcaMatrix *matrix);

int dca_matrix_multiply_classical(const DcaMatrix *left,
                                  const DcaMatrix *right,
                                  DcaMatrix *result,
                                  DcaOperationCounts *operations);

int dca_matrix_multiply_strassen(const DcaMatrix *left,
                                 const DcaMatrix *right,
                                 DcaMatrix *result,
                                 DcaOperationCounts *operations);

#endif
