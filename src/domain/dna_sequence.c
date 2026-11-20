#include "genome/domain/dna_sequence.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int dna_sequence_is_valid(const char *bases)
{
    size_t i;

    if (bases == NULL || bases[0] == '\0') {
        return 0;
    }

    for (i = 0; bases[i] != '\0'; ++i) {
        char base = bases[i];

        if (base != 'A' && base != 'C' && base != 'G' && base != 'T' &&
            base != 'a' && base != 'c' && base != 'g' && base != 't') {
            return 0;
        }
    }

    return 1;
}

int dna_sequence_create(DnaSequence *sequence,
                        const char *identifier,
                        const char *bases)
{
    size_t identifier_length;
    size_t bases_length;
    char *identifier_copy;
    char *bases_copy;

    if (sequence == NULL) {
        return 0;
    }

    if (sequence->identifier != NULL || sequence->bases != NULL ||
        sequence->length != 0) {
        return 0;
    }

    if (identifier == NULL || identifier[0] == '\0' ||
        !dna_sequence_is_valid(bases)) {
        return 0;
    }

    identifier_length = strlen(identifier);
    bases_length = strlen(bases);
    if (identifier_length == SIZE_MAX || bases_length == SIZE_MAX) {
        return 0;
    }

    identifier_copy = (char *)malloc(identifier_length + 1);
    if (identifier_copy == NULL) {
        return 0;
    }

    bases_copy = (char *)malloc(bases_length + 1);
    if (bases_copy == NULL) {
        free(identifier_copy);
        return 0;
    }

    memcpy(identifier_copy, identifier, identifier_length + 1);
    memcpy(bases_copy, bases, bases_length + 1);

    sequence->identifier = identifier_copy;
    sequence->bases = bases_copy;
    sequence->length = bases_length;

    return 1;
}

void dna_sequence_destroy(DnaSequence *sequence)
{
    if (sequence == NULL) {
        return;
    }

    free(sequence->identifier);
    free(sequence->bases);

    sequence->identifier = NULL;
    sequence->bases = NULL;
    sequence->length = 0;
}
