#include "genome/algorithms/divide_conquer.h"

#include <stdint.h>

ptrdiff_t dca_find_celebrity(const unsigned char *relationships,
                            size_t count,
                            DcaOperationCounts *operations)
{
    size_t candidate = 0;
    size_t person;

    if (operations != NULL) {
        operations->comparisons = 0;
        operations->additions = 0;
        operations->multiplications = 0;
    }
    if (count == 0 || relationships == NULL ||
        count > SIZE_MAX / count || count > (size_t)PTRDIFF_MAX) {
        return -1;
    }

    for (person = 1; person < count; ++person) {
        if (operations != NULL) {
            ++operations->comparisons;
        }
        if (relationships[candidate * count + person] != 0) {
            candidate = person;
        }
    }

    for (person = 0; person < count; ++person) {
        if (person == candidate) {
            continue;
        }
        if (operations != NULL) {
            ++operations->comparisons;
        }
        if (relationships[candidate * count + person] != 0) {
            return -1;
        }
        if (operations != NULL) {
            ++operations->comparisons;
        }
        if (relationships[person * count + candidate] == 0) {
            return -1;
        }
    }

    return (ptrdiff_t)candidate;
}
