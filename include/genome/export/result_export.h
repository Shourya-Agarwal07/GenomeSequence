#ifndef GENOME_EXPORT_RESULT_EXPORT_H
#define GENOME_EXPORT_RESULT_EXPORT_H

#include "genome/benchmark/benchmark.h"
#include "genome/domain/dna_sequence.h"

#include <stddef.h>

typedef struct {
    const char *algorithm;
    size_t result_count;
    int matches_reference;
} AlgorithmComparisonRecord;

int export_benchmark_records_csv(const char *path,
                                 const BenchmarkRecord *records,
                                 size_t count);
int export_benchmark_records_json(const char *path,
                                  const BenchmarkRecord *records,
                                  size_t count);
int export_sequence_statistics_csv(const char *path,
                                   const DnaSequence *sequences,
                                   size_t count);
int export_sequence_statistics_json(const char *path,
                                    const DnaSequence *sequences,
                                    size_t count);
int export_algorithm_comparisons_csv(
    const char *path,
    const AlgorithmComparisonRecord *records,
    size_t count);
int export_algorithm_comparisons_json(
    const char *path,
    const AlgorithmComparisonRecord *records,
    size_t count);

#endif
