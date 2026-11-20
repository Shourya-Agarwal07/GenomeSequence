#ifndef GENOME_ALGORITHMS_SEQUENCE_ALGORITHMS_H
#define GENOME_ALGORITHMS_SEQUENCE_ALGORITHMS_H

#include "genome/domain/dna_sequence.h"

#include <stddef.h>
#include <stdint.h>

typedef enum {
    SEQUENCE_SORT_BY_LENGTH,
    SEQUENCE_SORT_BY_IDENTIFIER
} SequenceSortCriterion;

typedef struct {
    size_t comparisons;
    size_t assignments;
    size_t swaps;
} SequenceAlgorithmStats;

typedef enum {
    SEQUENCE_SEARCH_ERROR = -1,
    SEQUENCE_SEARCH_NOT_FOUND = 0,
    SEQUENCE_SEARCH_FOUND = 1
} SequenceSearchResult;

int sequence_merge_sort(DnaSequence *records,
                        size_t count,
                        SequenceSortCriterion criterion,
                        SequenceAlgorithmStats *stats);
int sequence_quick_sort(DnaSequence *records,
                        size_t count,
                        SequenceSortCriterion criterion,
                        SequenceAlgorithmStats *stats);

SequenceSearchResult sequence_binary_search_iterative(
    const DnaSequence *records,
    size_t count,
    const DnaSequence *target,
    SequenceSortCriterion criterion,
    size_t *index,
    SequenceAlgorithmStats *stats);
SequenceSearchResult sequence_binary_search_recursive(
    const DnaSequence *records,
    size_t count,
    const DnaSequence *target,
    SequenceSortCriterion criterion,
    size_t *index,
    SequenceAlgorithmStats *stats);

#endif
