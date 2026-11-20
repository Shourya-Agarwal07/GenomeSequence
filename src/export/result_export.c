#include "genome/export/result_export.h"

#include "genome/sequence/sequence_statistics.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int finish_file(FILE *file, int successful)
{
    if (ferror(file)) {
        successful = 0;
    }
    if (fclose(file) != 0) {
        successful = 0;
    }
    return successful;
}

static int write_csv_field(FILE *file, const char *value)
{
    const unsigned char *cursor = (const unsigned char *)value;
    int quote = 0;

    while (*cursor != '\0') {
        if (*cursor == ',' || *cursor == '"' ||
            *cursor == '\r' || *cursor == '\n') {
            quote = 1;
            break;
        }
        ++cursor;
    }
    if (!quote) {
        return fputs(value, file) >= 0;
    }
    if (fputc('"', file) == EOF) {
        return 0;
    }
    cursor = (const unsigned char *)value;
    while (*cursor != '\0') {
        if (*cursor == '"' && fputc('"', file) == EOF) {
            return 0;
        }
        if (fputc(*cursor, file) == EOF) {
            return 0;
        }
        ++cursor;
    }
    return fputc('"', file) != EOF;
}

static int write_json_string(FILE *file, const char *value)
{
    const unsigned char *cursor = (const unsigned char *)value;

    if (fputc('"', file) == EOF) {
        return 0;
    }
    while (*cursor != '\0') {
        unsigned char byte = *cursor++;

        if (byte == '"' || byte == '\\') {
            if (fputc('\\', file) == EOF || fputc(byte, file) == EOF) {
                return 0;
            }
        } else if (byte == '\b') {
            if (fputs("\\b", file) < 0) {
                return 0;
            }
        } else if (byte == '\f') {
            if (fputs("\\f", file) < 0) {
                return 0;
            }
        } else if (byte == '\n') {
            if (fputs("\\n", file) < 0) {
                return 0;
            }
        } else if (byte == '\r') {
            if (fputs("\\r", file) < 0) {
                return 0;
            }
        } else if (byte == '\t') {
            if (fputs("\\t", file) < 0) {
                return 0;
            }
        } else if (byte < 0x20 || byte >= 0x80) {
            if (fprintf(file, "\\u%04X", (unsigned int)byte) < 0) {
                return 0;
            }
        } else if (fputc(byte, file) == EOF) {
            return 0;
        }
    }
    return fputc('"', file) != EOF;
}

static int benchmark_records_are_valid(const BenchmarkRecord *records,
                                       size_t count)
{
    size_t i;

    if (count != 0 && records == NULL) {
        return 0;
    }
    for (i = 0; i < count; ++i) {
        if (memchr(records[i].algorithm, '\0',
                   sizeof(records[i].algorithm)) == NULL ||
            memchr(records[i].dataset, '\0',
                   sizeof(records[i].dataset)) == NULL ||
            !isfinite(records[i].average_seconds) ||
            records[i].average_seconds < 0.0 ||
            records[i].repetitions == 0 ||
            (records[i].verified != 0 && records[i].verified != 1) ||
            (records[i].has_operation_count != 0 &&
             records[i].has_operation_count != 1) ||
            (records[i].has_measured_result_bytes != 0 &&
             records[i].has_measured_result_bytes != 1) ||
            (records[i].has_estimated_auxiliary_bytes != 0 &&
             records[i].has_estimated_auxiliary_bytes != 1) ||
            (records[i].has_result_count != 0 &&
             records[i].has_result_count != 1) ||
            (records[i].has_solution_value != 0 &&
             records[i].has_solution_value != 1)) {
            return 0;
        }
    }
    return 1;
}

static int sequence_records_are_valid(const DnaSequence *sequences,
                                      size_t count)
{
    size_t i;

    if (count != 0 && sequences == NULL) {
        return 0;
    }
    for (i = 0; i < count; ++i) {
        SequenceStatistics statistics;
        if (sequences[i].identifier == NULL ||
            (sequences[i].length != 0 && sequences[i].bases == NULL) ||
            !sequence_statistics_compute(sequences[i].bases,
                                         sequences[i].length,
                                         &statistics)) {
            return 0;
        }
    }
    return 1;
}

static int comparison_records_are_valid(
    const AlgorithmComparisonRecord *records,
    size_t count)
{
    size_t i;

    if (count != 0 && records == NULL) {
        return 0;
    }
    for (i = 0; i < count; ++i) {
        if (records[i].algorithm == NULL ||
            (records[i].matches_reference != 0 &&
             records[i].matches_reference != 1)) {
            return 0;
        }
    }
    return 1;
}

int export_benchmark_records_csv(const char *path,
                                 const BenchmarkRecord *records,
                                 size_t count)
{
    FILE *file;
    size_t i;
    int successful = 1;

    if (path == NULL || !benchmark_records_are_valid(records, count)) {
        return 0;
    }
    file = fopen(path, "wb");
    if (file == NULL) {
        return 0;
    }
    successful = fputs("algorithm,dataset,input_size,average_seconds,"
                       "repetitions,verified,operation_count,"
                       "measured_result_bytes,estimated_auxiliary_bytes,"
                       "result_count,solution_value\n", file) >= 0;
    for (i = 0; successful && i < count; ++i) {
        successful = write_csv_field(file, records[i].algorithm) &&
                     fputc(',', file) != EOF &&
                     write_csv_field(file, records[i].dataset) &&
                     fprintf(file, ",%zu,%.17g,%zu,%d,",
                             records[i].input_size,
                             records[i].average_seconds,
                             records[i].repetitions,
                             records[i].verified) >= 0;
        if (successful && records[i].has_operation_count) {
            successful = fprintf(file, "%llu",
                                 (unsigned long long)
                                     records[i].operation_count) >= 0;
        }
        if (successful) {
            successful = fputc(',', file) != EOF;
        }
        if (successful && records[i].has_measured_result_bytes) {
            successful = fprintf(file, "%zu",
                                 records[i].measured_result_bytes) >= 0;
        }
        if (successful) {
            successful = fputc(',', file) != EOF;
        }
        if (successful && records[i].has_estimated_auxiliary_bytes) {
            successful = fprintf(file, "%zu",
                                 records[i].estimated_auxiliary_bytes) >= 0;
        }
        if (successful) {
            successful = fputc(',', file) != EOF;
        }
        if (successful && records[i].has_result_count) {
            successful = fprintf(file, "%zu", records[i].result_count) >= 0;
        }
        if (successful) {
            successful = fputc(',', file) != EOF;
        }
        if (successful && records[i].has_solution_value) {
            successful = fprintf(file, "%lld",
                                 (long long)records[i].solution_value) >= 0;
        }
        if (successful) {
            successful = fputc('\n', file) != EOF;
        }
    }
    return finish_file(file, successful);
}

int export_benchmark_records_json(const char *path,
                                  const BenchmarkRecord *records,
                                  size_t count)
{
    FILE *file;
    size_t i;
    int successful;

    if (path == NULL || !benchmark_records_are_valid(records, count)) {
        return 0;
    }
    file = fopen(path, "wb");
    if (file == NULL) {
        return 0;
    }
    successful = fputs("{\"benchmarks\":[", file) >= 0;
    for (i = 0; successful && i < count; ++i) {
        if (i != 0) {
            successful = fputc(',', file) != EOF;
        }
        if (successful) {
            successful = fputs("{\"algorithm\":", file) >= 0 &&
                         write_json_string(file, records[i].algorithm) &&
                         fputs(",\"dataset\":", file) >= 0 &&
                         write_json_string(file, records[i].dataset) &&
                         fprintf(file, ",\"input_size\":%zu,"
                                 "\"average_seconds\":%.17g,"
                                 "\"repetitions\":%zu,\"verified\":%s,"
                                 "\"has_operation_count\":%s,"
                                 "\"operation_count\":%llu",
                                 records[i].input_size,
                                 records[i].average_seconds,
                                 records[i].repetitions,
                                 records[i].verified ? "true" : "false",
                                 records[i].has_operation_count
                                     ? "true"
                                     : "false",
                                 (unsigned long long)
                                     records[i].operation_count) >= 0;
            if (successful) {
                successful = fputs(",\"measured_result_bytes\":", file) >= 0;
                if (successful && records[i].has_measured_result_bytes) {
                    successful = fprintf(file, "%zu",
                                         records[i].measured_result_bytes) >= 0;
                } else if (successful) {
                    successful = fputs("null", file) >= 0;
                }
                successful = successful &&
                             fputs(",\"estimated_auxiliary_bytes\":", file) >= 0;
                if (successful && records[i].has_estimated_auxiliary_bytes) {
                    successful = fprintf(
                        file, "%zu", records[i].estimated_auxiliary_bytes) >= 0;
                } else if (successful) {
                    successful = fputs("null", file) >= 0;
                }
                successful = successful &&
                             fputs(",\"result_count\":", file) >= 0;
                if (successful && records[i].has_result_count) {
                    successful = fprintf(file, "%zu",
                                         records[i].result_count) >= 0;
                } else if (successful) {
                    successful = fputs("null", file) >= 0;
                }
                successful = successful &&
                             fputs(",\"solution_value\":", file) >= 0;
                if (successful && records[i].has_solution_value) {
                    successful = fprintf(file, "%lld",
                                         (long long)records[i].solution_value)
                                     >= 0;
                } else if (successful) {
                    successful = fputs("null", file) >= 0;
                }
                successful = successful && fputc('}', file) != EOF;
            }
        }
    }
    if (successful) {
        successful = fputs("]}\n", file) >= 0;
    }
    return finish_file(file, successful);
}

int export_sequence_statistics_csv(const char *path,
                                   const DnaSequence *sequences,
                                   size_t count)
{
    FILE *file;
    size_t i;
    int successful = 1;

    if (path == NULL || !sequence_records_are_valid(sequences, count)) {
        return 0;
    }
    file = fopen(path, "wb");
    if (file == NULL) {
        return 0;
    }
    successful = fputs("identifier,length,A,C,G,T,ambiguous,gc_percentage\n",
                       file) >= 0;
    for (i = 0; successful && i < count; ++i) {
        SequenceStatistics statistics;

        successful = sequence_statistics_compute(
                         sequences[i].bases, sequences[i].length,
                         &statistics) &&
                     write_csv_field(file, sequences[i].identifier) &&
                     fprintf(file, ",%zu,%zu,%zu,%zu,%zu,%zu,",
                             sequences[i].length,
                             statistics.a_count,
                             statistics.c_count,
                             statistics.g_count,
                             statistics.t_count,
                             statistics.ambiguous_count) >= 0;
        if (successful && statistics.gc_percentage_defined) {
            successful = fprintf(file, "%.17g",
                                 statistics.gc_percentage) >= 0;
        }
        if (successful) {
            successful = fputc('\n', file) != EOF;
        }
    }
    return finish_file(file, successful);
}

int export_sequence_statistics_json(const char *path,
                                    const DnaSequence *sequences,
                                    size_t count)
{
    FILE *file;
    size_t i;
    int successful;

    if (path == NULL || !sequence_records_are_valid(sequences, count)) {
        return 0;
    }
    file = fopen(path, "wb");
    if (file == NULL) {
        return 0;
    }
    successful = fputs("{\"sequences\":[", file) >= 0;
    for (i = 0; successful && i < count; ++i) {
        SequenceStatistics statistics;

        successful = sequence_statistics_compute(
            sequences[i].bases, sequences[i].length, &statistics);
        if (successful && i != 0) {
            successful = fputc(',', file) != EOF;
        }
        if (successful) {
            successful = fputs("{\"identifier\":", file) >= 0 &&
                         write_json_string(file, sequences[i].identifier) &&
                         fprintf(file,
                                 ",\"length\":%zu,\"A\":%zu,\"C\":%zu,"
                                 "\"G\":%zu,\"T\":%zu,"
                                 "\"ambiguous\":%zu,\"gc_percentage\":",
                                 sequences[i].length,
                                 statistics.a_count,
                                 statistics.c_count,
                                 statistics.g_count,
                                 statistics.t_count,
                                 statistics.ambiguous_count) >= 0;
        }
        if (successful && statistics.gc_percentage_defined) {
            successful = fprintf(file, "%.17g",
                                 statistics.gc_percentage) >= 0;
        } else if (successful) {
            successful = fputs("null", file) >= 0;
        }
        if (successful) {
            successful = fputc('}', file) != EOF;
        }
    }
    if (successful) {
        successful = fputs("]}\n", file) >= 0;
    }
    return finish_file(file, successful);
}

int export_algorithm_comparisons_csv(
    const char *path,
    const AlgorithmComparisonRecord *records,
    size_t count)
{
    FILE *file;
    size_t i;
    int successful = 1;

    if (path == NULL ||
        !comparison_records_are_valid(records, count)) {
        return 0;
    }
    file = fopen(path, "wb");
    if (file == NULL) {
        return 0;
    }
    successful = fputs("algorithm,result_count,matches_reference\n",
                       file) >= 0;
    for (i = 0; successful && i < count; ++i) {
        successful = write_csv_field(file, records[i].algorithm) &&
                     fprintf(file, ",%zu,%s\n",
                             records[i].result_count,
                             records[i].matches_reference ? "true" : "false")
                         >= 0;
    }
    return finish_file(file, successful);
}

int export_algorithm_comparisons_json(
    const char *path,
    const AlgorithmComparisonRecord *records,
    size_t count)
{
    FILE *file;
    size_t i;
    int successful;

    if (path == NULL ||
        !comparison_records_are_valid(records, count)) {
        return 0;
    }
    file = fopen(path, "wb");
    if (file == NULL) {
        return 0;
    }
    successful = fputs("{\"comparisons\":[", file) >= 0;
    for (i = 0; successful && i < count; ++i) {
        if (i != 0) {
            successful = fputc(',', file) != EOF;
        }
        if (successful) {
            successful = fputs("{\"algorithm\":", file) >= 0 &&
                         write_json_string(file, records[i].algorithm) &&
                         fprintf(file,
                                 ",\"result_count\":%zu,"
                                 "\"matches_reference\":%s}",
                                 records[i].result_count,
                                 records[i].matches_reference
                                     ? "true"
                                     : "false") >= 0;
        }
    }
    if (successful) {
        successful = fputs("]}\n", file) >= 0;
    }
    return finish_file(file, successful);
}
