#ifndef GENOME_STRING_MATCHING_STRING_MATCHING_H
#define GENOME_STRING_MATCHING_STRING_MATCHING_H

#include <stddef.h>

typedef struct {
    size_t *positions;
    size_t count;
    size_t capacity;
} StringMatchResult;

void string_match_result_init(StringMatchResult *result);
void string_match_result_destroy(StringMatchResult *result);

int string_match_naive(const char *text,
                       size_t text_length,
                       const char *pattern,
                       size_t pattern_length,
                       StringMatchResult *result);
int string_match_rabin_karp(const char *text,
                            size_t text_length,
                            const char *pattern,
                            size_t pattern_length,
                            StringMatchResult *result);
int string_match_kmp(const char *text,
                     size_t text_length,
                     const char *pattern,
                     size_t pattern_length,
                     StringMatchResult *result);
int string_match_finite_automaton(const char *text,
                                  size_t text_length,
                                  const char *pattern,
                                  size_t pattern_length,
                                  StringMatchResult *result);

#endif
