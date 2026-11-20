#ifndef GENOME_DOMAIN_DNA_SEQUENCE_H
#define GENOME_DOMAIN_DNA_SEQUENCE_H

#include <stddef.h>

typedef struct {
    char *identifier;
    char *bases;
    size_t length;
} DnaSequence;

int dna_sequence_is_valid(const char *bases);

int dna_sequence_create(DnaSequence *sequence,
                        const char *identifier,
                        const char *bases);

void dna_sequence_destroy(DnaSequence *sequence);

#endif
