#include "genome/algorithms/dna_alignment.h"

#include "genome/domain/dna_sequence.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

enum {
    TRACE_STOP = 0,
    TRACE_DIAGONAL = 1,
    TRACE_UP = 2,
    TRACE_LEFT = 3
};

static int checked_add_int64(int64_t left, int64_t right, int64_t *sum)
{
    if ((right > 0 && left > INT64_MAX - right) ||
        (right < 0 && left < INT64_MIN - right)) {
        return 0;
    }
    *sum = left + right;
    return 1;
}

static char normalized_base(char base)
{
    if (base >= 'a' && base <= 'z') {
        return (char)(base - 'a' + 'A');
    }
    return base;
}

static int allocation_dimensions(size_t rows,
                                size_t columns,
                                size_t *cell_count)
{
    if (columns != 0 && rows > SIZE_MAX / columns) {
        return 0;
    }
    *cell_count = rows * columns;
    return *cell_count <= SIZE_MAX / sizeof(int64_t);
}

void dna_alignment_result_destroy(DnaAlignmentResult *result)
{
    if (result == NULL) {
        return;
    }
    free(result->aligned_first);
    free(result->aligned_second);
    free(result->match_line);
    memset(result, 0, sizeof(*result));
}

int dna_alignment_compute(const char *first,
                          const char *second,
                          DnaAlignmentType type,
                          int64_t match_score,
                          int64_t mismatch_penalty,
                          int64_t gap_penalty,
                          DnaAlignmentResult *result)
{
    size_t first_length;
    size_t second_length;
    size_t rows;
    size_t columns;
    size_t cell_count;
    size_t max_alignment_length;
    size_t alignment_capacity;
    int64_t *scores = NULL;
    unsigned char *trace = NULL;
    DnaAlignmentResult computed = {0};
    size_t best_row = 0;
    size_t best_column = 0;
    int64_t best_score = 0;
    size_t i;
    size_t j;
    size_t aligned_length = 0;
    int success = 0;

    if (first == NULL || second == NULL || result == NULL ||
        !dna_sequence_is_valid(first) || !dna_sequence_is_valid(second) ||
        (type != DNA_ALIGNMENT_GLOBAL && type != DNA_ALIGNMENT_LOCAL)) {
        return 0;
    }

    first_length = strlen(first);
    second_length = strlen(second);
    if (first_length == SIZE_MAX || second_length == SIZE_MAX ||
        first_length > SIZE_MAX - second_length) {
        return 0;
    }
    max_alignment_length = first_length + second_length;
    if (max_alignment_length == SIZE_MAX) {
        return 0;
    }
    alignment_capacity = max_alignment_length + 1;
    rows = first_length + 1;
    columns = second_length + 1;
    if (!allocation_dimensions(rows, columns, &cell_count) ||
        cell_count > SIZE_MAX / sizeof(*trace)) {
        return 0;
    }

    scores = (int64_t *)malloc(cell_count * sizeof(*scores));
    trace = (unsigned char *)malloc(cell_count * sizeof(*trace));
    computed.aligned_first = (char *)malloc(alignment_capacity);
    computed.aligned_second = (char *)malloc(alignment_capacity);
    computed.match_line = (char *)malloc(alignment_capacity);
    if (scores == NULL || trace == NULL || computed.aligned_first == NULL ||
        computed.aligned_second == NULL || computed.match_line == NULL) {
        goto cleanup;
    }

    memset(trace, TRACE_STOP, cell_count * sizeof(*trace));
    if (type == DNA_ALIGNMENT_GLOBAL) {
        scores[0] = 0;
        for (i = 1; i < rows; ++i) {
            if (!checked_add_int64(scores[(i - 1) * columns], gap_penalty,
                                   &scores[i * columns])) {
                goto cleanup;
            }
            trace[i * columns] = TRACE_UP;
        }
        for (j = 1; j < columns; ++j) {
            if (!checked_add_int64(scores[j - 1], gap_penalty, &scores[j])) {
                goto cleanup;
            }
            trace[j] = TRACE_LEFT;
        }
    } else {
        memset(scores, 0, cell_count * sizeof(*scores));
    }

    for (i = 1; i < rows; ++i) {
        for (j = 1; j < columns; ++j) {
            size_t cell = i * columns + j;
            int64_t diagonal;
            int64_t up;
            int64_t left;
            int64_t best;
            unsigned char direction;
            int64_t substitution =
                normalized_base(first[i - 1]) ==
                        normalized_base(second[j - 1])
                    ? match_score
                    : mismatch_penalty;

            if (!checked_add_int64(scores[(i - 1) * columns + j - 1],
                                   substitution, &diagonal) ||
                !checked_add_int64(scores[(i - 1) * columns + j],
                                   gap_penalty, &up) ||
                !checked_add_int64(scores[i * columns + j - 1],
                                   gap_penalty, &left)) {
                goto cleanup;
            }

            best = diagonal;
            direction = TRACE_DIAGONAL;
            if (up > best) {
                best = up;
                direction = TRACE_UP;
            }
            if (left > best) {
                best = left;
                direction = TRACE_LEFT;
            }
            if (type == DNA_ALIGNMENT_LOCAL && best <= 0) {
                best = 0;
                direction = TRACE_STOP;
            }
            scores[cell] = best;
            trace[cell] = direction;
            if (type == DNA_ALIGNMENT_LOCAL && best > best_score) {
                best_score = best;
                best_row = i;
                best_column = j;
            }
        }
    }

    if (type == DNA_ALIGNMENT_GLOBAL) {
        best_row = first_length;
        best_column = second_length;
        best_score = scores[best_row * columns + best_column];
    }
    computed.type = type;
    computed.score = best_score;
    computed.first_end = best_row;
    computed.second_end = best_column;

    i = best_row;
    j = best_column;
    while ((type == DNA_ALIGNMENT_GLOBAL && (i != 0 || j != 0)) ||
           (type == DNA_ALIGNMENT_LOCAL && i != 0 && j != 0 &&
            trace[i * columns + j] != TRACE_STOP)) {
        unsigned char direction = trace[i * columns + j];

        if (direction == TRACE_DIAGONAL) {
            char first_base = normalized_base(first[i - 1]);
            char second_base = normalized_base(second[j - 1]);
            computed.aligned_first[aligned_length] = first_base;
            computed.aligned_second[aligned_length] = second_base;
            computed.match_line[aligned_length] =
                first_base == second_base ? '|' : ' ';
            if (first_base == second_base) {
                ++computed.matches;
            } else {
                ++computed.mismatches;
            }
            --i;
            --j;
        } else if (direction == TRACE_UP && i != 0) {
            computed.aligned_first[aligned_length] =
                normalized_base(first[i - 1]);
            computed.aligned_second[aligned_length] = '-';
            computed.match_line[aligned_length] = ' ';
            ++computed.gaps;
            --i;
        } else if (direction == TRACE_LEFT && j != 0) {
            computed.aligned_first[aligned_length] = '-';
            computed.aligned_second[aligned_length] =
                normalized_base(second[j - 1]);
            computed.match_line[aligned_length] = ' ';
            ++computed.gaps;
            --j;
        } else {
            goto cleanup;
        }
        ++aligned_length;
    }
    computed.first_start = i;
    computed.second_start = j;
    computed.length = aligned_length;
    computed.aligned_first[aligned_length] = '\0';
    computed.aligned_second[aligned_length] = '\0';
    computed.match_line[aligned_length] = '\0';
    for (i = 0; i < aligned_length / 2; ++i) {
        char temporary;
        size_t opposite = aligned_length - i - 1;

        temporary = computed.aligned_first[i];
        computed.aligned_first[i] = computed.aligned_first[opposite];
        computed.aligned_first[opposite] = temporary;
        temporary = computed.aligned_second[i];
        computed.aligned_second[i] = computed.aligned_second[opposite];
        computed.aligned_second[opposite] = temporary;
        temporary = computed.match_line[i];
        computed.match_line[i] = computed.match_line[opposite];
        computed.match_line[opposite] = temporary;
    }

    dna_alignment_result_destroy(result);
    *result = computed;
    memset(&computed, 0, sizeof(computed));
    success = 1;

cleanup:
    free(scores);
    free(trace);
    dna_alignment_result_destroy(&computed);
    return success;
}
