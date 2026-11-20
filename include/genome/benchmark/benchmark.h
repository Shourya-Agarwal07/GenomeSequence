#ifndef GENOME_BENCHMARK_BENCHMARK_H
#define GENOME_BENCHMARK_BENCHMARK_H

#include <stddef.h>
#include <stdint.h>

#define BENCHMARK_ALGORITHM_NAME_CAPACITY 64
#define BENCHMARK_DATASET_NAME_CAPACITY 32
#define BENCHMARK_SUITE_MAX_RECORDS 42

typedef struct {
    char algorithm[BENCHMARK_ALGORITHM_NAME_CAPACITY];
    char dataset[BENCHMARK_DATASET_NAME_CAPACITY];
    size_t input_size;
    double average_seconds;
    size_t repetitions;
    int verified;
    uint64_t operation_count;
    int has_operation_count;
    size_t measured_result_bytes;
    int has_measured_result_bytes;
    size_t estimated_auxiliary_bytes;
    int has_estimated_auxiliary_bytes;
    size_t result_count;
    int has_result_count;
    int64_t solution_value;
    int has_solution_value;
} BenchmarkRecord;

typedef int (*BenchmarkCallback)(void *context);
typedef void (*BenchmarkCleanup)(void *context);
typedef int (*BenchmarkOperationCount)(void *context, uint64_t *count);

int benchmark_measure(const char *algorithm,
                      size_t input_size,
                      size_t repetitions,
                      BenchmarkCallback prepare,
                      BenchmarkCallback operation,
                      BenchmarkCallback verify,
                      BenchmarkCleanup cleanup,
                      BenchmarkOperationCount operation_count,
                      void *context,
                      BenchmarkRecord *record);

int benchmark_run_suite(BenchmarkRecord *records,
                        size_t record_capacity,
                        size_t *record_count);

#endif
