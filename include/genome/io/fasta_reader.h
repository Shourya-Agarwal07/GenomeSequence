#ifndef GENOME_IO_FASTA_READER_H
#define GENOME_IO_FASTA_READER_H

#include <stddef.h>
#include "genome/domain/dna_sequence.h"

int fasta_read_file(const char *path,
                    DnaSequence **sequences_out,
                    size_t *count_out);

void fasta_sequences_destroy(DnaSequence *sequences, size_t count);

#endif
