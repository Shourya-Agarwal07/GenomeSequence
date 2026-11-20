#include "genome/io/fasta_reader.h"

#include <ctype.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int reserve_buffer(char **buffer, size_t *capacity, size_t required)
{
    size_t new_capacity;
    char *grown;

    if (required <= *capacity)
    {
        return 1;
    }

    new_capacity = *capacity == 0 ? 128 : *capacity;
    while (new_capacity < required)
    {
        if (new_capacity > SIZE_MAX / 2)
        {
            new_capacity = required;
            break;
        }
        new_capacity *= 2;
    }

    grown = (char *)realloc(*buffer, new_capacity);
    if (grown == NULL)
    {
        return 0;
    }

    *buffer = grown;
    *capacity = new_capacity;
    return 1;
}

static int read_line(FILE *file,
                     char **line,
                     size_t *line_length,
                     size_t *line_capacity)
{
    int character;

    *line_length = 0;
    while ((character = fgetc(file)) != EOF)
    {
        if (character == '\n')
        {
            break;
        }

        if (*line_length > SIZE_MAX - 2 ||
            !reserve_buffer(line, line_capacity, *line_length + 2))
        {
            return -1;
        }

        (*line)[(*line_length)++] = (char)character;
    }

    if (character == EOF && ferror(file))
    {
        return -1;
    }
    if (character == EOF && *line_length == 0)
    {
        return 0;
    }

    if (*line_length > 0 && (*line)[*line_length - 1] == '\r')
    {
        --(*line_length);
    }
    (*line)[*line_length] = '\0';
    return 1;
}

static int is_whitespace_only(const char *text, size_t length)
{
    size_t i;

    for (i = 0; i < length; ++i)
    {
        if (!isspace((unsigned char)text[i]))
        {
            return 0;
        }
    }
    return 1;
}

static int append_bases(char **bases,
                        size_t *length,
                        size_t *capacity,
                        const char *line,
                        size_t line_length)
{
    size_t required;

    if (*length == SIZE_MAX || line_length > SIZE_MAX - *length)
    {
        return 0;
    }
    required = *length + line_length;
    if (required == SIZE_MAX)
    {
        return 0;
    }
    ++required;

    if (!reserve_buffer(bases, capacity, required))
    {
        return 0;
    }

    memcpy(*bases + *length, line, line_length);
    *length += line_length;
    (*bases)[*length] = '\0';
    return 1;
}

static int append_sequence(DnaSequence **sequences,
                           size_t *count,
                           size_t *capacity,
                           const char *identifier,
                           const char *bases)
{
    DnaSequence item = {NULL, NULL, 0};

    if (!dna_sequence_create(&item, identifier, bases))
    {
        return 0;
    }

    if (*count == *capacity)
    {
        size_t new_capacity;
        DnaSequence *grown;

        if (*capacity == 0)
        {
            new_capacity = 4;
        }
        else
        {
            if (*capacity > SIZE_MAX / 2)
            {
                dna_sequence_destroy(&item);
                return 0;
            }
            new_capacity = *capacity * 2;
        }
        if (new_capacity > SIZE_MAX / sizeof(**sequences))
        {
            dna_sequence_destroy(&item);
            return 0;
        }

        grown = (DnaSequence *)realloc(
            *sequences, new_capacity * sizeof(**sequences));
        if (grown == NULL)
        {
            dna_sequence_destroy(&item);
            return 0;
        }

        *sequences = grown;
        *capacity = new_capacity;
    }

    (*sequences)[(*count)++] = item;
    return 1;
}

void fasta_sequences_destroy(DnaSequence *sequences, size_t count)
{
    size_t i;

    if (sequences == NULL)
    {
        return;
    }

    for (i = 0; i < count; ++i)
    {
        dna_sequence_destroy(&sequences[i]);
    }

    free(sequences);
}

int fasta_read_file(const char *path,
                    DnaSequence **sequences_out,
                    size_t *count_out)
{
    FILE *file;
    char *line = NULL;
    size_t line_length = 0;
    size_t line_capacity = 0;
    char *current_identifier = NULL;
    char *current_bases = NULL;
    size_t bases_length = 0;
    size_t bases_capacity = 0;
    DnaSequence *sequences = NULL;
    size_t count = 0;
    size_t capacity = 0;
    int success = 0;
    int line_status;

    if (sequences_out != NULL)
    {
        *sequences_out = NULL;
    }
    if (count_out != NULL)
    {
        *count_out = 0;
    }
    if (path == NULL || sequences_out == NULL || count_out == NULL)
    {
        return 0;
    }

    file = fopen(path, "rb");
    if (file == NULL)
    {
        return 0;
    }

    while ((line_status = read_line(file, &line, &line_length,
                                    &line_capacity)) == 1)
    {
        if (line_length == 0)
        {
            continue;
        }
        if (memchr(line, '\0', line_length) != NULL)
        {
            goto cleanup;
        }

        if (line[0] == '>')
        {
            if (current_identifier != NULL)
            {
                if (bases_length == 0 ||
                    !append_sequence(&sequences, &count, &capacity,
                                     current_identifier, current_bases))
                {
                    goto cleanup;
                }

                bases_length = 0;
                if (current_bases != NULL)
                {
                    current_bases[0] = '\0';
                }
            }

            free(current_identifier);
            current_identifier = NULL;
            if (line_length == 1 ||
                is_whitespace_only(line + 1, line_length - 1))
            {
                goto cleanup;
            }
            current_identifier = (char *)malloc(line_length);
            if (current_identifier == NULL)
            {
                goto cleanup;
            }
            memcpy(current_identifier, line + 1, line_length - 1);
            current_identifier[line_length - 1] = '\0';
        }
        else
        {
            if (current_identifier == NULL)
            {
                goto cleanup;
            }

            if (!append_bases(&current_bases, &bases_length,
                              &bases_capacity, line, line_length))
            {
                goto cleanup;
            }
        }
    }

    if (line_status < 0 || ferror(file))
    {
        goto cleanup;
    }

    if (current_identifier != NULL)
    {
        if (bases_length == 0 ||
            !append_sequence(&sequences, &count, &capacity,
                             current_identifier, current_bases))
        {
            goto cleanup;
        }
    }

    if (count == 0)
    {
        goto cleanup;
    }

    *sequences_out = sequences;
    *count_out = count;
    sequences = NULL;
    count = 0;
    success = 1;

cleanup:
    fclose(file);
    free(line);
    free(current_identifier);
    free(current_bases);
    fasta_sequences_destroy(sequences, count);
    return success;
}
