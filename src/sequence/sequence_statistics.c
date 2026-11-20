#include "genome/sequence/sequence_statistics.h"

int sequence_statistics_compute(const char *bases,
                                size_t length,
                                SequenceStatistics *statistics)
{
    SequenceStatistics computed = {0, 0, 0, 0, 0, 0.0, 0};
    size_t i;

    if (statistics == NULL) {
        return 0;
    }
    *statistics = computed;
    if (length != 0 && bases == NULL) {
        return 0;
    }

    for (i = 0; i < length; ++i) {
        switch (bases[i]) {
        case 'A':
        case 'a':
            ++computed.a_count;
            break;
        case 'C':
        case 'c':
            ++computed.c_count;
            break;
        case 'G':
        case 'g':
            ++computed.g_count;
            break;
        case 'T':
        case 't':
            ++computed.t_count;
            break;
        default:
            ++computed.ambiguous_count;
            break;
        }
    }

    if (length != 0 && computed.ambiguous_count == 0) {
        computed.gc_percentage =
            100.0 * (double)(computed.g_count + computed.c_count) /
            (double)length;
        computed.gc_percentage_defined = 1;
    }
    *statistics = computed;
    return 1;
}
