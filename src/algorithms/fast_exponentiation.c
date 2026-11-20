#include "genome/algorithms/divide_conquer.h"

#include <math.h>
#include <stdint.h>

int dca_fast_power(double base,
                   int exponent,
                   double *result,
                   DcaOperationCounts *operations)
{
    uint64_t remaining;
    double factor = base;
    double value = 1.0;

    if (operations != NULL) {
        operations->comparisons = 0;
        operations->additions = 0;
        operations->multiplications = 0;
    }
    if (result != NULL) {
        *result = 0.0;
    }
    if (result == NULL || !isfinite(base)) {
        return 0;
    }
    if (exponent < 0) {
        if (base == 0.0) {
            return 0;
        }
        factor = 1.0 / base;
    }

    if (exponent < 0) {
        remaining = (uint64_t)(-(exponent + 1));
        ++remaining;
    } else {
        remaining = (uint64_t)exponent;
    }
    while (remaining != 0) {
        if ((remaining & UINT64_C(1)) != 0) {
            value *= factor;
            if (operations != NULL) {
                ++operations->multiplications;
            }
            if (!isfinite(value)) {
                return 0;
            }
        }

        remaining >>= 1;
        if (remaining != 0) {
            factor *= factor;
            if (operations != NULL) {
                ++operations->multiplications;
            }
            if (!isfinite(factor)) {
                return 0;
            }
        }
    }

    *result = value;
    return 1;
}
