#include "genome/string_matching/string_matching.h"

#include <stdint.h>
#include <stdlib.h>
#include <string.h>

static int result_append(StringMatchResult *result, size_t position)
{
    if (result->count == result->capacity) {
        size_t new_capacity;
        size_t *new_positions;

        if (result->count == SIZE_MAX) {
            return 0;
        }
        if (result->capacity == 0) {
            new_capacity = 8;
        } else if (result->capacity <= SIZE_MAX / 2) {
            new_capacity = result->capacity * 2;
        } else {
            new_capacity = SIZE_MAX;
        }
        if (new_capacity < result->count + 1) {
            new_capacity = result->count + 1;
        }
        if (new_capacity > SIZE_MAX / sizeof(*new_positions)) {
            return 0;
        }
        new_positions = (size_t *)realloc(
            result->positions, new_capacity * sizeof(*new_positions));
        if (new_positions == NULL) {
            return 0;
        }
        result->positions = new_positions;
        result->capacity = new_capacity;
    }

    result->positions[result->count] = position;
    ++result->count;
    return 1;
}

static int append_empty_pattern_matches(size_t text_length,
                                        StringMatchResult *result)
{
    size_t position = 0;

    if (text_length == SIZE_MAX) {
        return 0;
    }
    for (;;) {
        if (!result_append(result, position)) {
            return 0;
        }
        if (position == text_length) {
            return 1;
        }
        ++position;
    }
}

static void replace_result(StringMatchResult *result,
                           StringMatchResult *replacement)
{
    string_match_result_destroy(result);
    *result = *replacement;
    string_match_result_init(replacement);
}

void string_match_result_init(StringMatchResult *result)
{
    if (result != NULL) {
        result->positions = NULL;
        result->count = 0;
        result->capacity = 0;
    }
}

void string_match_result_destroy(StringMatchResult *result)
{
    if (result == NULL) {
        return;
    }
    free(result->positions);
    string_match_result_init(result);
}

int string_match_naive(const char *text,
                       size_t text_length,
                       const char *pattern,
                       size_t pattern_length,
                       StringMatchResult *result)
{
    StringMatchResult matches = {NULL, 0, 0};
    size_t window_count;
    size_t start;

    if (text == NULL || pattern == NULL || result == NULL) {
        return 0;
    }
    if (pattern_length == 0) {
        if (!append_empty_pattern_matches(text_length, &matches)) {
            string_match_result_destroy(&matches);
            return 0;
        }
        replace_result(result, &matches);
        return 1;
    }
    if (pattern_length > text_length) {
        replace_result(result, &matches);
        return 1;
    }

    window_count = text_length - pattern_length + 1;
    for (start = 0; start < window_count; ++start) {
        size_t offset = 0;

        while (offset < pattern_length &&
               text[start + offset] == pattern[offset]) {
            ++offset;
        }
        if (offset == pattern_length &&
            !result_append(&matches, start)) {
            string_match_result_destroy(&matches);
            return 0;
        }
    }

    replace_result(result, &matches);
    return 1;
}

int string_match_rabin_karp(const char *text,
                            size_t text_length,
                            const char *pattern,
                            size_t pattern_length,
                            StringMatchResult *result)
{
    const uint64_t base = UINT64_C(256);
    const uint64_t modulus = UINT64_C(1000000007);
    StringMatchResult matches = {NULL, 0, 0};
    uint64_t pattern_hash = 0;
    uint64_t window_hash = 0;
    uint64_t highest_place = 1;
    size_t window_count;
    size_t i;
    size_t start;

    if (text == NULL || pattern == NULL || result == NULL) {
        return 0;
    }
    if (pattern_length == 0) {
        if (!append_empty_pattern_matches(text_length, &matches)) {
            string_match_result_destroy(&matches);
            return 0;
        }
        replace_result(result, &matches);
        return 1;
    }
    if (pattern_length > text_length) {
        replace_result(result, &matches);
        return 1;
    }

    window_count = text_length - pattern_length + 1;
    for (i = 0; i < pattern_length; ++i) {
        pattern_hash =
            (pattern_hash * base + (unsigned char)pattern[i]) % modulus;
        window_hash =
            (window_hash * base + (unsigned char)text[i]) % modulus;
        if (i + 1 < pattern_length) {
            highest_place = (highest_place * base) % modulus;
        }
    }

    for (start = 0; start < window_count; ++start) {
        if (pattern_hash == window_hash &&
            memcmp(text + start, pattern, pattern_length) == 0 &&
            !result_append(&matches, start)) {
            string_match_result_destroy(&matches);
            return 0;
        }

        if (start + 1 < window_count) {
            uint64_t outgoing =
                ((uint64_t)(unsigned char)text[start] * highest_place) %
                modulus;
            uint64_t incoming =
                (unsigned char)text[start + pattern_length];

            window_hash = (window_hash + modulus - outgoing) % modulus;
            window_hash = (window_hash * base + incoming) % modulus;
        }
    }

    replace_result(result, &matches);
    return 1;
}

int string_match_kmp(const char *text,
                     size_t text_length,
                     const char *pattern,
                     size_t pattern_length,
                     StringMatchResult *result)
{
    StringMatchResult matches = {NULL, 0, 0};
    size_t *lps = NULL;
    size_t prefix_length = 0;
    size_t i;
    int success = 0;

    if (text == NULL || pattern == NULL || result == NULL) {
        return 0;
    }
    if (pattern_length == 0) {
        if (!append_empty_pattern_matches(text_length, &matches)) {
            string_match_result_destroy(&matches);
            return 0;
        }
        replace_result(result, &matches);
        return 1;
    }
    if (pattern_length > text_length) {
        replace_result(result, &matches);
        return 1;
    }
    if (pattern_length > SIZE_MAX / sizeof(*lps)) {
        return 0;
    }
    lps = (size_t *)malloc(pattern_length * sizeof(*lps));
    if (lps == NULL) {
        return 0;
    }

    lps[0] = 0;
    i = 1;
    while (i < pattern_length) {
        if (pattern[i] == pattern[prefix_length]) {
            ++prefix_length;
            lps[i] = prefix_length;
            ++i;
        } else if (prefix_length != 0) {
            prefix_length = lps[prefix_length - 1];
        } else {
            lps[i] = 0;
            ++i;
        }
    }

    {
        size_t text_index = 0;
        size_t pattern_index = 0;

        while (text_index < text_length) {
            if (text[text_index] == pattern[pattern_index]) {
                ++text_index;
                ++pattern_index;
                if (pattern_index == pattern_length) {
                    if (!result_append(&matches,
                                       text_index - pattern_length)) {
                        goto cleanup;
                    }
                    pattern_index = lps[pattern_index - 1];
                }
            } else if (pattern_index != 0) {
                pattern_index = lps[pattern_index - 1];
            } else {
                ++text_index;
            }
        }
    }

    replace_result(result, &matches);
    success = 1;

cleanup:
    free(lps);
    if (!success) {
        string_match_result_destroy(&matches);
    }
    return success;
}

int string_match_finite_automaton(const char *text,
                                  size_t text_length,
                                  const char *pattern,
                                  size_t pattern_length,
                                  StringMatchResult *result)
{
    const size_t alphabet_size = 256;
    StringMatchResult matches = {NULL, 0, 0};
    size_t *lps = NULL;
    size_t *transitions = NULL;
    size_t transition_states;
    size_t transition_count;
    size_t prefix_length = 0;
    size_t i;
    size_t state;
    int success = 0;

    if (text == NULL || pattern == NULL || result == NULL) {
        return 0;
    }
    if (pattern_length == 0) {
        if (!append_empty_pattern_matches(text_length, &matches)) {
            string_match_result_destroy(&matches);
            return 0;
        }
        replace_result(result, &matches);
        return 1;
    }
    if (pattern_length > text_length) {
        replace_result(result, &matches);
        return 1;
    }
    if (pattern_length == SIZE_MAX ||
        pattern_length > SIZE_MAX / sizeof(*lps)) {
        return 0;
    }
    transition_states = pattern_length + 1;
    if (transition_states > SIZE_MAX / alphabet_size) {
        return 0;
    }
    transition_count = transition_states * alphabet_size;
    if (transition_count > SIZE_MAX / sizeof(*transitions)) {
        return 0;
    }
    lps = (size_t *)malloc(pattern_length * sizeof(*lps));
    transitions = (size_t *)malloc(transition_count *
                                   sizeof(*transitions));
    if (lps == NULL || transitions == NULL) {
        goto cleanup;
    }

    lps[0] = 0;
    i = 1;
    while (i < pattern_length) {
        if ((unsigned char)pattern[i] ==
            (unsigned char)pattern[prefix_length]) {
            ++prefix_length;
            lps[i] = prefix_length;
            ++i;
        } else if (prefix_length != 0) {
            prefix_length = lps[prefix_length - 1];
        } else {
            lps[i] = 0;
            ++i;
        }
    }

    for (i = 0; i < alphabet_size; ++i) {
        transitions[i] =
            (unsigned char)pattern[0] == (unsigned char)i ? 1 : 0;
    }
    for (state = 1; state < pattern_length; ++state) {
        size_t fallback = lps[state - 1];
        size_t row = state * alphabet_size;
        size_t fallback_row = fallback * alphabet_size;

        for (i = 0; i < alphabet_size; ++i) {
            transitions[row + i] =
                (unsigned char)pattern[state] == (unsigned char)i
                    ? state + 1
                    : transitions[fallback_row + i];
        }
    }
    {
        size_t row = pattern_length * alphabet_size;
        size_t fallback_row = lps[pattern_length - 1] * alphabet_size;

        for (i = 0; i < alphabet_size; ++i) {
            transitions[row + i] = transitions[fallback_row + i];
        }
    }

    state = 0;
    for (i = 0; i < text_length; ++i) {
        state = transitions[state * alphabet_size +
                            (unsigned char)text[i]];
        if (state == pattern_length &&
            !result_append(&matches, i + 1 - pattern_length)) {
            goto cleanup;
        }
    }

    replace_result(result, &matches);
    success = 1;

cleanup:
    free(transitions);
    free(lps);
    if (!success) {
        string_match_result_destroy(&matches);
    }
    return success;
}
