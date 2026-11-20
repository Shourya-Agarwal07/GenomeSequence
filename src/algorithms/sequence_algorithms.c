#include "genome/algorithms/sequence_algorithms.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int criterion_is_valid(SequenceSortCriterion criterion)
{
    return criterion == SEQUENCE_SORT_BY_LENGTH ||
           criterion == SEQUENCE_SORT_BY_IDENTIFIER;
}

static int records_are_valid(const DnaSequence *records,
                             size_t count,
                             SequenceSortCriterion criterion)
{
    size_t i;

    if (!criterion_is_valid(criterion) || (count != 0 && records == NULL)) {
        return 0;
    }

    if (criterion == SEQUENCE_SORT_BY_IDENTIFIER) {
        for (i = 0; i < count; ++i) {
            if (records[i].identifier == NULL) {
                return 0;
            }
        }
    }

    return 1;
}

static int compare_records(const DnaSequence *left,
                           const DnaSequence *right,
                           SequenceSortCriterion criterion)
{
    if (criterion == SEQUENCE_SORT_BY_LENGTH) {
        if (left->length < right->length) {
            return -1;
        }
        if (left->length > right->length) {
            return 1;
        }
        return 0;
    }

    return strcmp(left->identifier, right->identifier);
}

static void reset_stats(SequenceAlgorithmStats *stats)
{
    if (stats != NULL) {
        stats->comparisons = 0;
        stats->assignments = 0;
        stats->swaps = 0;
    }
}

static void merge_ranges(DnaSequence *records,
                         DnaSequence *temporary,
                         size_t start,
                         size_t middle,
                         size_t end,
                         SequenceSortCriterion criterion,
                         SequenceAlgorithmStats *stats)
{
    size_t left = start;
    size_t right = middle;
    size_t output = start;
    size_t i;

    while (left < middle && right < end) {
        if (stats != NULL) {
            ++stats->comparisons;
        }
        if (compare_records(&records[left], &records[right], criterion) <= 0) {
            temporary[output++] = records[left++];
        } else {
            temporary[output++] = records[right++];
        }
        if (stats != NULL) {
            ++stats->assignments;
        }
    }

    while (left < middle) {
        temporary[output++] = records[left++];
        if (stats != NULL) {
            ++stats->assignments;
        }
    }
    while (right < end) {
        temporary[output++] = records[right++];
        if (stats != NULL) {
            ++stats->assignments;
        }
    }

    for (i = start; i < end; ++i) {
        records[i] = temporary[i];
        if (stats != NULL) {
            ++stats->assignments;
        }
    }
}

int sequence_merge_sort(DnaSequence *records,
                        size_t count,
                        SequenceSortCriterion criterion,
                        SequenceAlgorithmStats *stats)
{
    DnaSequence *temporary;
    size_t width;

    reset_stats(stats);
    if (!records_are_valid(records, count, criterion)) {
        return 0;
    }
    if (count < 2) {
        return 1;
    }
    if (count > SIZE_MAX / sizeof(*temporary)) {
        return 0;
    }

    temporary = (DnaSequence *)malloc(count * sizeof(*temporary));
    if (temporary == NULL) {
        return 0;
    }

    for (width = 1; width < count;) {
        size_t start;

        for (start = 0; start < count;) {
            size_t middle = start + (width < count - start
                                         ? width
                                         : count - start);
            size_t end = middle + (width < count - middle
                                       ? width
                                       : count - middle);

            if (middle < end) {
                merge_ranges(records, temporary, start, middle, end,
                             criterion, stats);
            }
            start = end;
        }

        if (width > count / 2) {
            break;
        }
        width *= 2;
    }

    free(temporary);
    return 1;
}

static void swap_records(DnaSequence *left,
                         DnaSequence *right,
                         SequenceAlgorithmStats *stats)
{
    DnaSequence temporary;

    if (left == right) {
        return;
    }

    temporary = *left;
    *left = *right;
    *right = temporary;
    if (stats != NULL) {
        ++stats->swaps;
        stats->assignments += 3;
    }
}

static void quick_sort_range(DnaSequence *records,
                             size_t start,
                             size_t end,
                             SequenceSortCriterion criterion,
                             SequenceAlgorithmStats *stats)
{
    while (end - start > 1) {
        size_t pivot_index = start + (end - start) / 2;
        DnaSequence pivot = records[pivot_index];
        size_t less = start;
        size_t current = start;
        size_t greater = end;

        if (stats != NULL) {
            ++stats->assignments;
        }

        while (current < greater) {
            int comparison;

            if (stats != NULL) {
                ++stats->comparisons;
            }
            comparison = compare_records(&records[current], &pivot, criterion);
            if (comparison < 0) {
                swap_records(&records[current], &records[less], stats);
                ++current;
                ++less;
            } else if (comparison > 0) {
                --greater;
                swap_records(&records[current], &records[greater], stats);
            } else {
                ++current;
            }
        }

        if (less - start < end - greater) {
            quick_sort_range(records, start, less, criterion, stats);
            start = greater;
        } else {
            quick_sort_range(records, greater, end, criterion, stats);
            end = less;
        }
    }
}

int sequence_quick_sort(DnaSequence *records,
                        size_t count,
                        SequenceSortCriterion criterion,
                        SequenceAlgorithmStats *stats)
{
    reset_stats(stats);
    if (!records_are_valid(records, count, criterion)) {
        return 0;
    }

    if (count > 1) {
        quick_sort_range(records, 0, count, criterion, stats);
    }
    return 1;
}

static int search_arguments_are_valid(const DnaSequence *records,
                                      size_t count,
                                      const DnaSequence *target,
                                      SequenceSortCriterion criterion)
{
    return target != NULL &&
           records_are_valid(records, count, criterion) &&
           (criterion != SEQUENCE_SORT_BY_IDENTIFIER ||
            target->identifier != NULL);
}

static SequenceSearchResult binary_search_recursive_range(
    const DnaSequence *records,
    size_t start,
    size_t end,
    const DnaSequence *target,
    SequenceSortCriterion criterion,
    size_t *index,
    SequenceAlgorithmStats *stats)
{
    size_t middle;
    int comparison;

    if (start >= end) {
        return SEQUENCE_SEARCH_NOT_FOUND;
    }

    middle = start + (end - start) / 2;
    if (stats != NULL) {
        ++stats->comparisons;
    }
    comparison = compare_records(&records[middle], target, criterion);
    if (comparison == 0) {
        *index = middle;
        return SEQUENCE_SEARCH_FOUND;
    }
    if (comparison < 0) {
        return binary_search_recursive_range(records, middle + 1, end,
                                             target, criterion, index, stats);
    }
    return binary_search_recursive_range(records, start, middle, target,
                                         criterion, index, stats);
}

SequenceSearchResult sequence_binary_search_iterative(
    const DnaSequence *records,
    size_t count,
    const DnaSequence *target,
    SequenceSortCriterion criterion,
    size_t *index,
    SequenceAlgorithmStats *stats)
{
    size_t start = 0;
    size_t end = count;

    reset_stats(stats);
    if (index != NULL) {
        *index = SIZE_MAX;
    }
    if (index == NULL ||
        !search_arguments_are_valid(records, count, target, criterion)) {
        return SEQUENCE_SEARCH_ERROR;
    }

    while (start < end) {
        size_t middle = start + (end - start) / 2;
        int comparison;

        if (stats != NULL) {
            ++stats->comparisons;
        }
        comparison = compare_records(&records[middle], target, criterion);
        if (comparison == 0) {
            *index = middle;
            return SEQUENCE_SEARCH_FOUND;
        }
        if (comparison < 0) {
            start = middle + 1;
        } else {
            end = middle;
        }
    }

    return SEQUENCE_SEARCH_NOT_FOUND;
}

SequenceSearchResult sequence_binary_search_recursive(
    const DnaSequence *records,
    size_t count,
    const DnaSequence *target,
    SequenceSortCriterion criterion,
    size_t *index,
    SequenceAlgorithmStats *stats)
{
    reset_stats(stats);
    if (index != NULL) {
        *index = SIZE_MAX;
    }
    if (index == NULL ||
        !search_arguments_are_valid(records, count, target, criterion)) {
        return SEQUENCE_SEARCH_ERROR;
    }

    return binary_search_recursive_range(records, 0, count, target, criterion,
                                         index, stats);
}
