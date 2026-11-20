#ifndef GENOME_ALGORITHMS_DNA_ALIGNMENT_H
#define GENOME_ALGORITHMS_DNA_ALIGNMENT_H

#include <stddef.h>
#include <stdint.h>

typedef enum {
    DNA_ALIGNMENT_GLOBAL,
    DNA_ALIGNMENT_LOCAL
} DnaAlignmentType;

typedef struct {
    DnaAlignmentType type;
    int64_t score;
    char *aligned_first;
    char *aligned_second;
    char *match_line;
    size_t length;
    size_t matches;
    size_t mismatches;
    size_t gaps;
    size_t first_start;
    size_t first_end;
    size_t second_start;
    size_t second_end;
} DnaAlignmentResult;

int dna_alignment_compute(const char *first,
                          const char *second,
                          DnaAlignmentType type,
                          int64_t match_score,
                          int64_t mismatch_penalty,
                          int64_t gap_penalty,
                          DnaAlignmentResult *result);

void dna_alignment_result_destroy(DnaAlignmentResult *result);

#endif
