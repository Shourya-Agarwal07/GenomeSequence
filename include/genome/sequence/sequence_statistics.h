#ifndef GENOME_SEQUENCE_SEQUENCE_STATISTICS_H
#define GENOME_SEQUENCE_SEQUENCE_STATISTICS_H

#include <stddef.h>

typedef struct {
    size_t a_count;
    size_t c_count;
    size_t g_count;
    size_t t_count;
    size_t ambiguous_count;
    double gc_percentage;
    int gc_percentage_defined;
} SequenceStatistics;

int sequence_statistics_compute(const char *bases,
                                size_t length,
                                SequenceStatistics *statistics);

#endif
