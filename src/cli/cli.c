#include "genome/cli/cli.h"

#include "genome/algorithms/backtracking_algorithms.h"
#include "genome/algorithms/dna_alignment.h"
#include "genome/algorithms/divide_conquer.h"
#include "genome/algorithms/dynamic_programming.h"
#include "genome/algorithms/greedy_algorithms.h"
#include "genome/algorithms/sequence_algorithms.h"
#include "genome/benchmark/benchmark.h"
#include "genome/domain/dna_sequence.h"
#include "genome/graph/bipartite_matching.h"
#include "genome/graph/ford_fulkerson.h"
#include "genome/io/fasta_reader.h"
#include "genome/sequence/sequence_statistics.h"
#include "genome/string_matching/string_matching.h"

#include <errno.h>
#include <inttypes.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef _WIN32
#include <io.h>
#include <windows.h>
#else
#include <unistd.h>
#endif

#define CLI_LINE_CAPACITY 2048
#define CLI_ALIGNMENT_HISTORY_CAPACITY 10
#define STRESS_TEXT_LENGTH 1000000
#define STRESS_SEQUENCE_COUNT 512
#define STRESS_SEQUENCE_LENGTH 128

typedef struct {
    DnaAlignmentType type;
    int64_t score;
    size_t matches;
    size_t mismatches;
    size_t gaps;
    char first[CLI_LINE_CAPACITY];
    char second[CLI_LINE_CAPACITY];
} CliAlignmentHistoryEntry;

typedef struct {
    int terminal_enabled;
    int color_enabled;
    int64_t match_score;
    int64_t mismatch_penalty;
    int64_t gap_penalty;
    CliAlignmentHistoryEntry alignment_history[CLI_ALIGNMENT_HISTORY_CAPACITY];
    size_t alignment_history_count;
#ifdef _WIN32
    HANDLE console_handle;
    DWORD original_console_mode;
    int console_mode_changed;
#endif
} CliContext;

static int output_is_terminal(FILE *output)
{
    if (output != stdout) {
        return 0;
    }
#ifdef _WIN32
    return _isatty(_fileno(output)) != 0;
#else
    return isatty(fileno(output)) != 0;
#endif
}

static void configure_terminal(CliContext *context, FILE *output)
{
    const char *no_color = getenv("NO_COLOR");
    const char *term = getenv("TERM");

    context->terminal_enabled = output_is_terminal(output);
    if (!context->terminal_enabled) {
        return;
    }
#ifdef _WIN32
    context->console_handle = GetStdHandle(STD_OUTPUT_HANDLE);
    if (context->console_handle != INVALID_HANDLE_VALUE &&
        GetConsoleMode(context->console_handle,
                       &context->original_console_mode)) {
        DWORD mode = context->original_console_mode |
                     ENABLE_VIRTUAL_TERMINAL_PROCESSING;
        if (SetConsoleMode(context->console_handle, mode)) {
            context->console_mode_changed = 1;
            context->color_enabled = 1;
        } else {
            context->terminal_enabled = 0;
        }
    }
#else
    context->color_enabled = 1;
#endif
    if (no_color != NULL) {
        context->color_enabled = 0;
    }
    if (term != NULL && strcmp(term, "dumb") == 0) {
        context->terminal_enabled = 0;
        context->color_enabled = 0;
    }
}

static void restore_terminal(CliContext *context)
{
#ifdef _WIN32
    if (context->console_mode_changed) {
        (void)SetConsoleMode(context->console_handle,
                             context->original_console_mode);
        context->console_mode_changed = 0;
    }
#else
    (void)context;
#endif
}

static void clear_screen(const CliContext *context, FILE *output)
{
    if (context->terminal_enabled) {
        fputs("\033[2J\033[H", output);
    }
}

static void print_heading(FILE *output,
                          const CliContext *context,
                          const char *heading)
{
    fputs("\n============================================================\n", output);
    if (context->color_enabled) {
        fputs("\033[1;36m", output);
    }
    fprintf(output, "%s\n", heading);
    if (context->color_enabled) {
        fputs("\033[0m", output);
    }
    fputs("============================================================\n", output);
}

static void print_error(FILE *output,
                        const CliContext *context,
                        const char *message)
{
    if (context->color_enabled) {
        fputs("\033[1;31m", output);
    }
    fprintf(output, "Error: %s\n", message);
    if (context->color_enabled) {
        fputs("\033[0m", output);
    }
}

static void print_success(FILE *output,
                          const CliContext *context,
                          const char *message)
{
    if (context->color_enabled) {
        fputs("\033[1;32m", output);
    }
    fprintf(output, "%s\n", message);
    if (context->color_enabled) {
        fputs("\033[0m", output);
    }
}

static int read_line(FILE *input,
                     FILE *output,
                     const char *prompt,
                     char *buffer,
                     size_t capacity)
{
    char *newline;

    if (prompt != NULL) {
        fputs(prompt, output);
        fflush(output);
    }
    if (fgets(buffer, (int)capacity, input) == NULL) {
        return 0;
    }
    newline = strchr(buffer, '\n');
    if (newline == NULL && !feof(input)) {
        int character;
        while ((character = fgetc(input)) != '\n' && character != EOF) {
        }
        fputs("Input is too long; please try again.\n", output);
        return -1;
    }
    if (newline != NULL) {
        *newline = '\0';
    }
    newline = strchr(buffer, '\r');
    if (newline != NULL) {
        *newline = '\0';
    }
    return 1;
}

static int parse_int64(const char *text, int64_t *value)
{
    char *end;
    intmax_t parsed;

    errno = 0;
    parsed = strtoimax(text, &end, 10);
    if (text == end || errno == ERANGE || parsed < INT64_MIN ||
        parsed > INT64_MAX) {
        return 0;
    }
    while (*end == ' ' || *end == '\t') {
        ++end;
    }
    if (*end != '\0') {
        return 0;
    }
    *value = (int64_t)parsed;
    return 1;
}

static int parse_size(const char *text, size_t *value)
{
    char *end;
    uintmax_t parsed;

    if (*text == '-') {
        return 0;
    }
    errno = 0;
    parsed = strtoumax(text, &end, 10);
    if (text == end || errno == ERANGE || parsed > SIZE_MAX) {
        return 0;
    }
    while (*end == ' ' || *end == '\t') {
        ++end;
    }
    if (*end != '\0') {
        return 0;
    }
    *value = (size_t)parsed;
    return 1;
}

static int read_menu_choice(FILE *input,
                            FILE *output,
                            size_t maximum_choice,
                            size_t *choice)
{
    char line[64];
    int status = read_line(input, output, "Enter your choice: ",
                           line, sizeof(line));

    if (status != 1) {
        return status;
    }
    if (!parse_size(line, choice) || *choice > maximum_choice) {
        return -1;
    }
    return 1;
}

static int read_size_value(FILE *input,
                           FILE *output,
                           const char *prompt,
                           size_t *value)
{
    char line[CLI_LINE_CAPACITY];
    int status = read_line(input, output, prompt, line, sizeof(line));

    if (status != 1) {
        return status;
    }
    if (!parse_size(line, value)) {
        fputs("Please enter a nonnegative whole number.\n", output);
        return -1;
    }
    return 1;
}

static int prompt_fasta(FILE *input,
                        FILE *output,
                        DnaSequence **sequences,
                        size_t *count)
{
    char path[CLI_LINE_CAPACITY];
    int status = read_line(input, output,
                           "FASTA path (blank uses data/samples/demo.fasta): ",
                           path, sizeof(path));

    if (status != 1) {
        return 0;
    }
    if (path[0] == '\0') {
        (void)snprintf(path, sizeof(path), "data/samples/demo.fasta");
    }
    if (!fasta_read_file(path, sequences, count)) {
        fprintf(output, "Could not load valid FASTA data from '%s'.\n", path);
        return 0;
    }
    fprintf(output, "Loaded %zu sequence(s).\n", *count);
    return 1;
}

static void print_statistics(FILE *output,
                             const char *identifier,
                             const char *bases,
                             size_t length)
{
    SequenceStatistics statistics;

    if (!sequence_statistics_compute(bases, length, &statistics)) {
        fputs("Could not calculate sequence statistics.\n", output);
        return;
    }
    fprintf(output,
            "%s: length=%zu A=%zu C=%zu G=%zu T=%zu ambiguous=%zu",
            identifier, length, statistics.a_count, statistics.c_count,
            statistics.g_count, statistics.t_count,
            statistics.ambiguous_count);
    if (statistics.gc_percentage_defined) {
        fprintf(output, " GC=%.2f%%\n", statistics.gc_percentage);
    } else {
        fputs(" GC=undefined\n", output);
    }
}

static void print_string_matches(FILE *output,
                                 const char *name,
                                 const StringMatchResult *result)
{
    size_t i;

    fprintf(output, "%s: %zu occurrence(s)", name, result->count);
    for (i = 0; i < result->count; ++i) {
        fprintf(output, "%s%zu", i == 0 ? " at " : ", ",
                result->positions[i]);
    }
    fputc('\n', output);
}

static int results_are_equal(const StringMatchResult *left,
                             const StringMatchResult *right)
{
    return left->count == right->count &&
           (left->count == 0 ||
            memcmp(left->positions, right->positions,
                   left->count * sizeof(*left->positions)) == 0);
}

static int run_dna_statistics(FILE *input, FILE *output)
{
    char bases[CLI_LINE_CAPACITY];
    DnaSequence sequence = {NULL, NULL, 0};
    int status;

    status = read_line(input, output, "DNA bases: ", bases, sizeof(bases));
    if (status != 1) {
        return status == 0;
    }
    if (!dna_sequence_create(&sequence, "interactive", bases)) {
        fputs("Invalid DNA sequence; use only A, C, G, and T.\n", output);
        return 1;
    }
    print_statistics(output, sequence.identifier, sequence.bases,
                     sequence.length);
    dna_sequence_destroy(&sequence);
    return 1;
}

static int run_fasta_display(FILE *input, FILE *output)
{
    DnaSequence *sequences = NULL;
    size_t count = 0;
    size_t i;

    if (!prompt_fasta(input, output, &sequences, &count)) {
        return 1;
    }
    for (i = 0; i < count; ++i) {
        fprintf(output, "%zu. %s | length=%zu | %s\n", i + 1,
                sequences[i].identifier, sequences[i].length,
                sequences[i].bases);
        print_statistics(output, sequences[i].identifier, sequences[i].bases,
                         sequences[i].length);
    }
    fasta_sequences_destroy(sequences, count);
    return 1;
}

static int run_sort_and_search(FILE *input, FILE *output)
{
    DnaSequence *sequences = NULL;
    DnaSequence target = {NULL, NULL, 0};
    size_t count = 0;
    size_t key_choice;
    size_t search_choice;
    size_t query_length;
    size_t index = SIZE_MAX;
    SequenceSortCriterion criterion;
    SequenceSearchResult result;
    char query[CLI_LINE_CAPACITY];
    size_t i;

    if (!prompt_fasta(input, output, &sequences, &count)) {
        return 1;
    }
    if (read_size_value(input, output, "Sort key: 1=length, 2=identifier: ",
                        &key_choice) != 1 ||
        (key_choice != 1 && key_choice != 2)) {
        fputs("Invalid sort key.\n", output);
        fasta_sequences_destroy(sequences, count);
        return 1;
    }
    criterion = key_choice == 1 ? SEQUENCE_SORT_BY_LENGTH
                                : SEQUENCE_SORT_BY_IDENTIFIER;
    if (!sequence_quick_sort(sequences, count, criterion, NULL)) {
        fputs("Could not sort the loaded records.\n", output);
        fasta_sequences_destroy(sequences, count);
        return 1;
    }
    for (i = 0; i < count; ++i) {
        fprintf(output, "%zu. %s | length=%zu\n", i + 1,
                sequences[i].identifier, sequences[i].length);
    }

    if (criterion == SEQUENCE_SORT_BY_LENGTH) {
        if (read_size_value(input, output, "Length to search for: ",
                            &query_length) != 1) {
            fasta_sequences_destroy(sequences, count);
            return 1;
        }
        target.length = query_length;
    } else {
        if (read_line(input, output, "Identifier to search for: ", query,
                      sizeof(query)) != 1) {
            fasta_sequences_destroy(sequences, count);
            return 1;
        }
        target.identifier = query;
    }
    if (read_size_value(input, output,
                        "Search method: 1=iterative, 2=recursive: ",
                        &search_choice) != 1 ||
        (search_choice != 1 && search_choice != 2)) {
        fputs("Invalid search method.\n", output);
        fasta_sequences_destroy(sequences, count);
        return 1;
    }
    result = search_choice == 1
                 ? sequence_binary_search_iterative(
                       sequences, count, &target, criterion, &index, NULL)
                 : sequence_binary_search_recursive(
                       sequences, count, &target, criterion, &index, NULL);
    if (result == SEQUENCE_SEARCH_FOUND) {
        fprintf(output, "Found %s at sorted index %zu.\n",
                sequences[index].identifier, index);
    } else if (result == SEQUENCE_SEARCH_NOT_FOUND) {
        fputs("No matching record was found.\n", output);
    } else {
        fputs("Search rejected the provided data.\n", output);
    }
    fasta_sequences_destroy(sequences, count);
    return 1;
}

static int run_divide_and_conquer(FILE *output)
{
    const int64_t values[] = {-2, 1, -3, 4, -1, 2, 1, -5, 4};
    MaximumSubarrayResult result;

    if (!dca_maximum_subarray(values,
                              sizeof(values) / sizeof(values[0]),
                              &result, NULL)) {
        fputs("Maximum-subarray example failed.\n", output);
        return 1;
    }
    fprintf(output,
            "Maximum subarray: sum=%" PRId64 " indices=[%zu,%zu]. "
            "Divide and conquer: O(n log n) time, O(log n) stack.\n",
            result.sum, result.start_index, result.end_index);
    return 1;
}

static int run_greedy(FILE *output)
{
    const Activity activities[] = {
        {1, 4, 1}, {3, 5, 2}, {0, 6, 3}, {5, 7, 4}, {8, 9, 5}
    };
    Activity selected[5];
    size_t selected_count = 0;
    size_t i;

    if (!greedy_activity_selection(activities, 5, selected, 5,
                                   &selected_count)) {
        fputs("Activity-selection example failed.\n", output);
        return 1;
    }
    fprintf(output, "Selected %zu compatible activities:", selected_count);
    for (i = 0; i < selected_count; ++i) {
        fprintf(output, " %zu", selected[i].identifier);
    }
    fputs("\nGreedy activity selection: O(n log n) time.\n", output);
    return 1;
}

static int run_lcs_workflow(FILE *input, FILE *output)
{
    char first[CLI_LINE_CAPACITY];
    char second[CLI_LINE_CAPACITY];
    size_t lcs_length = 0;
    char *lcs = NULL;

    if (read_line(input, output, "First DNA sequence: ", first,
                  sizeof(first)) != 1 ||
        read_line(input, output, "Second DNA sequence: ", second,
                  sizeof(second)) != 1) {
        return 1;
    }
    if (!dna_sequence_is_valid(first) || !dna_sequence_is_valid(second)) {
        fputs("Both inputs must contain only A, C, G, and T.\n", output);
        return 1;
    }
    if (!dp_longest_common_subsequence(first, second, &lcs_length, &lcs)) {
        fputs("LCS computation failed.\n", output);
        return 1;
    }
    fprintf(output,
            "Bytewise LCS length=%zu sequence=%s. This is a sequence "
            "comparison demonstration, not biological alignment.\n",
            lcs_length, lcs);
    free(lcs);
    return 1;
}

static int read_alignment_score(FILE *input,
                                FILE *output,
                                const CliContext *context,
                                const char *description,
                                int64_t *value)
{
    char prompt[192];
    char line[CLI_LINE_CAPACITY];
    int status;

    for (;;) {
        (void)snprintf(prompt, sizeof(prompt), "%s [%+" PRId64 "]: ",
                       description, *value);
        status = read_line(input, output, prompt, line, sizeof(line));
        if (status == 0) {
            return 0;
        }
        if (status < 0) {
            continue;
        }
        if (line[0] == '\0') {
            return 1;
        }
        if (parse_int64(line, value)) {
            return 1;
        }
        print_error(output, context,
                    "Enter a signed integer, or press Enter to keep the "
                    "displayed value.");
    }
}

static void print_alignment_scoring(FILE *output, const CliContext *context)
{
    print_heading(output, context, "CURRENT ALIGNMENT SCORING PARAMETERS");
    fprintf(output,
            "Match reward     : %+" PRId64 "\n"
            "Mismatch penalty : %" PRId64 "\n"
            "Gap penalty      : %" PRId64 "\n",
            context->match_score, context->mismatch_penalty,
            context->gap_penalty);
    fputs("These values are used as the defaults for the next alignment.\n",
          output);
}

static void print_alignment_history(FILE *output, const CliContext *context)
{
    size_t i;

    print_heading(output, context, "ALIGNMENT HISTORY (THIS SESSION)");
    if (context->alignment_history_count == 0) {
        fputs("No alignments have been calculated in this session.\n",
              output);
        return;
    }
    for (i = 0; i < context->alignment_history_count; ++i) {
        const CliAlignmentHistoryEntry *entry =
            &context->alignment_history[i];
        const char *algorithm =
            entry->type == DNA_ALIGNMENT_GLOBAL
                ? "Needleman-Wunsch"
                : "Smith-Waterman";

        fprintf(output,
                "[%02zu] %-16s score=%" PRId64
                " matches=%zu mismatches=%zu gaps=%zu\n"
                "     %s  vs  %s\n",
                i + 1, algorithm, entry->score, entry->matches,
                entry->mismatches, entry->gaps, entry->first,
                entry->second);
    }
}

static void remember_alignment(CliContext *context,
                               const char *first,
                               const char *second,
                               const DnaAlignmentResult *result)
{
    CliAlignmentHistoryEntry *entry;

    if (context->alignment_history_count <
        CLI_ALIGNMENT_HISTORY_CAPACITY) {
        entry = &context->alignment_history[
            context->alignment_history_count++];
    } else {
        memmove(context->alignment_history,
                context->alignment_history + 1,
                (CLI_ALIGNMENT_HISTORY_CAPACITY - 1) *
                    sizeof(context->alignment_history[0]));
        entry = &context->alignment_history[
            CLI_ALIGNMENT_HISTORY_CAPACITY - 1];
    }
    entry->type = result->type;
    entry->score = result->score;
    entry->matches = result->matches;
    entry->mismatches = result->mismatches;
    entry->gaps = result->gaps;
    (void)snprintf(entry->first, sizeof(entry->first), "%s", first);
    (void)snprintf(entry->second, sizeof(entry->second), "%s", second);
}

static void run_alignment_calculation(FILE *input,
                                      FILE *output,
                                      CliContext *context,
                                      DnaAlignmentType type)
{
    char first[CLI_LINE_CAPACITY];
    char second[CLI_LINE_CAPACITY];
    DnaAlignmentResult result = {0};
    const char *algorithm = type == DNA_ALIGNMENT_GLOBAL
                                ? "Needleman-Wunsch"
                                : "Smith-Waterman";
    const char *alignment_type = type == DNA_ALIGNMENT_GLOBAL
                                     ? "Global Alignment"
                                     : "Local Alignment";
    double identity;

    if (read_line(input, output,
                  "Enter first DNA sequence (A, C, G, T): ",
                  first, sizeof(first)) != 1 ||
        read_line(input, output,
                  "Enter second DNA sequence (A, C, G, T): ",
                  second, sizeof(second)) != 1) {
        return;
    }
    if (!dna_sequence_is_valid(first) || !dna_sequence_is_valid(second)) {
        print_error(output, context,
                    "Both sequences must be non-empty and contain only A, "
                    "C, G, and T.");
        return;
    }
    if (!read_alignment_score(
            input, output, context,
            "Enter match reward (points for identical DNA characters)",
            &context->match_score) ||
        !read_alignment_score(
            input, output, context,
            "Enter mismatch penalty (points deducted for different "
            "characters)",
            &context->mismatch_penalty) ||
        !read_alignment_score(
            input, output, context,
            "Enter gap penalty (points deducted for inserting a gap)",
            &context->gap_penalty)) {
        return;
    }
    if (!dna_alignment_compute(first, second, type, context->match_score,
                               context->mismatch_penalty,
                               context->gap_penalty, &result)) {
        print_error(output, context,
                    "Alignment failed because of score overflow or "
                    "insufficient memory.");
        return;
    }

    remember_alignment(context, first, second, &result);
    identity = result.length == 0
                   ? 0.0
                   : 100.0 * (double)result.matches / (double)result.length;
    print_heading(output, context, "ALIGNMENT RESULT");
    fprintf(output,
            "Algorithm       : %s\n"
            "Alignment Type  : %s\n"
            "Sequence 1      : %s\n"
            "Sequence 2      : %s\n\n"
            "Aligned Sequence 1 : %s\n"
            "                     %s\n"
            "Aligned Sequence 2 : %s\n\n",
            algorithm, alignment_type, first, second,
            result.aligned_first, result.match_line,
            result.aligned_second);
    fprintf(output,
            "Alignment Score : %" PRId64 "\n"
            "Matches         : %zu\n"
            "Mismatches      : %zu\n"
            "Gaps            : %zu\n"
            "Identity        : %.2f%%\n",
            result.score,
            result.matches, result.mismatches, result.gaps, identity);
    if (type == DNA_ALIGNMENT_LOCAL && result.length != 0) {
        fprintf(output,
                "Matching region in sequence 1: %zu-%zu (1-based, inclusive)\n"
                "Matching region in sequence 2: %zu-%zu (1-based, inclusive)\n",
                result.first_start + 1, result.first_end,
                result.second_start + 1, result.second_end);
    } else if (type == DNA_ALIGNMENT_LOCAL) {
        fputs("No positive-scoring local region was found.\n", output);
    }
    fputs("============================================================\n", output);
    print_success(output, context, "Alignment calculated successfully.");
    dna_alignment_result_destroy(&result);
}

static int run_dna_alignment(FILE *input,
                             FILE *output,
                             CliContext *context)
{
    int first_display = 1;

    for (;;) {
        size_t choice;
        int status;

        if (first_display) {
            clear_screen(context, output);
            first_display = 0;
        }
        print_heading(output, context, "DNA SEQUENCE ALIGNMENT");
        fputs("Compare and align DNA sequences using dynamic programming.\n\n"
              "[01] Needleman-Wunsch\n"
              "     Global alignment of complete DNA sequences.\n"
              "[02] Smith-Waterman\n"
              "     Local alignment of matching DNA regions.\n"
              "[03] View Scoring Parameters\n"
              "[04] Alignment History\n"
              "[00] Return to Main Menu\n\n",
              output);
        status = read_menu_choice(input, output, 4, &choice);
        if (status == 0) {
            return 1;
        }
        if (status < 0) {
            print_error(output, context,
                        "Invalid alignment choice. Choose 01, 02, 03, 04, "
                        "or 00.");
            continue;
        }
        if (choice == 0) {
            return 1;
        }
        if (choice == 1) {
            run_alignment_calculation(input, output, context,
                                      DNA_ALIGNMENT_GLOBAL);
        } else if (choice == 2) {
            run_alignment_calculation(input, output, context,
                                      DNA_ALIGNMENT_LOCAL);
        } else if (choice == 3) {
            print_alignment_scoring(output, context);
        } else {
            print_alignment_history(output, context);
        }
    }
}

static int run_branch_and_bound(FILE *input, FILE *output)
{
    const size_t weights[] = {2, 3, 4, 5, 9};
    const int64_t values[] = {3, 4, 5, 8, 10};
    const char *names[] = {"DP", "Backtracking", "FIFO B&B", "LIFO B&B",
                           "Least-cost B&B"};
    KnapsackSearchResult results[4] = {{0, NULL, 0}, {0, NULL, 0},
                                       {0, NULL, 0}, {0, NULL, 0}};
    size_t capacity;
    int64_t dp_value = 0;
    int (*solvers[])(const size_t *, const int64_t *, size_t, size_t,
                     KnapsackSearchResult *) = {
        backtracking_knapsack_01,
        branch_bound_fifo_knapsack_01,
        branch_bound_lifo_knapsack_01,
        branch_bound_least_cost_knapsack_01
    };
    size_t i;

    if (read_size_value(input, output, "Knapsack capacity (0-23): ",
                        &capacity) != 1) {
        return 1;
    }
    if (capacity > 23 ||
        !dp_knapsack_01(weights, values, 5, capacity, &dp_value)) {
        fputs("Invalid capacity.\n", output);
        return 1;
    }
    fprintf(output, "%s result: %" PRId64 "\n", names[0], dp_value);
    for (i = 0; i < 4; ++i) {
        if (!solvers[i](weights, values, 5, capacity, &results[i])) {
            fputs("A search implementation failed.\n", output);
            while (i != 0) {
                knapsack_search_result_destroy(&results[--i]);
            }
            return 1;
        }
        fprintf(output, "%s result: %" PRId64 "\n", names[i + 1],
                results[i].maximum_value);
    }
    for (i = 0; i < 4; ++i) {
        knapsack_search_result_destroy(&results[i]);
    }
    return 1;
}

static int run_graph_flow(FILE *output)
{
    const FordFulkersonEdge edges[] = {
        {0, 1, 16}, {0, 2, 13}, {1, 2, 10}, {2, 1, 4},
        {1, 3, 12}, {3, 2, 9}, {2, 4, 14}, {4, 3, 7},
        {3, 5, 20}, {4, 5, 4}
    };
    FordFulkersonNetwork network = {6, 10, edges};
    FordFulkersonResult result = {0};
    const BipartiteEdge matching_edges[] = {
        {0, 0}, {0, 1}, {1, 0}, {2, 1}, {2, 2}
    };
    BipartiteGraph graph = {3, 3, 5, matching_edges};
    BipartiteMatchingResult matching = {0};

    if (!ford_fulkerson_max_flow(&network, 0, 5, &result)) {
        fputs("Flow example failed.\n", output);
        return 1;
    }
    fprintf(output, "Textbook flow network maximum flow: %" PRId64
                    " (cut capacity %" PRId64 ").\n",
            result.max_flow, result.cut_capacity);
    ford_fulkerson_result_destroy(&result);
    if (!bipartite_maximum_matching(&graph, &matching)) {
        fputs("Bipartite matching example failed.\n", output);
        return 1;
    }
    fprintf(output, "Maximum bipartite matching size: %zu.\n",
            matching.matching_size);
    bipartite_matching_result_destroy(&matching);
    return 1;
}

static int run_string_search(FILE *input, FILE *output)
{
    DnaSequence *sequences = NULL;
    size_t count = 0;
    size_t sequence_number;
    size_t algorithm_choice;
    char pattern[CLI_LINE_CAPACITY];
    size_t pattern_length;
    StringMatchResult results[4];
    const char *names[] = {"Naive", "Rabin-Karp", "KMP",
                           "Finite Automaton"};
    size_t i;
    int compare_all;

    if (!prompt_fasta(input, output, &sequences, &count)) {
        return 1;
    }
    if (read_size_value(input, output, "Sequence number: ",
                        &sequence_number) != 1 ||
        sequence_number == 0 || sequence_number > count) {
        fputs("Invalid sequence number.\n", output);
        fasta_sequences_destroy(sequences, count);
        return 1;
    }
    if (read_line(input, output, "DNA pattern: ", pattern,
                  sizeof(pattern)) != 1 ||
        !dna_sequence_is_valid(pattern)) {
        fputs("Pattern must be a nonempty DNA sequence.\n", output);
        fasta_sequences_destroy(sequences, count);
        return 1;
    }
    pattern_length = strlen(pattern);
    if (read_size_value(input, output,
                        "Algorithm: 1=Naive, 2=Rabin-Karp, 3=KMP, "
                        "4=compare, 5=Finite Automaton: ",
                        &algorithm_choice) != 1 ||
        algorithm_choice < 1 || algorithm_choice > 5) {
        fputs("Invalid string-matching algorithm.\n", output);
        fasta_sequences_destroy(sequences, count);
        return 1;
    }

    for (i = 0; i < 4; ++i) {
        string_match_result_init(&results[i]);
    }
    compare_all = algorithm_choice == 4;
    if (compare_all || algorithm_choice == 1) {
        if (!string_match_naive(sequences[sequence_number - 1].bases,
                                sequences[sequence_number - 1].length,
                                pattern, pattern_length, &results[0])) {
            goto failure;
        }
    }
    if (compare_all || algorithm_choice == 2) {
        if (!string_match_rabin_karp(
                sequences[sequence_number - 1].bases,
                sequences[sequence_number - 1].length, pattern,
                pattern_length, &results[1])) {
            goto failure;
        }
    }
    if (compare_all || algorithm_choice == 3) {
        if (!string_match_kmp(sequences[sequence_number - 1].bases,
                              sequences[sequence_number - 1].length,
                              pattern, pattern_length, &results[2])) {
            goto failure;
        }
    }
    if (compare_all || algorithm_choice == 5) {
        if (!string_match_finite_automaton(
                sequences[sequence_number - 1].bases,
                sequences[sequence_number - 1].length, pattern,
                pattern_length, &results[3])) {
            goto failure;
        }
    }
    if (compare_all) {
        int equal = results_are_equal(&results[0], &results[1]) &&
                    results_are_equal(&results[0], &results[2]) &&
                    results_are_equal(&results[0], &results[3]);
        fprintf(output, "All algorithms agree: %s.\n",
                equal ? "yes" : "no");
        for (i = 0; i < 4; ++i) {
            print_string_matches(output, names[i], &results[i]);
        }
    } else {
        size_t selected = algorithm_choice == 5 ? 3 : algorithm_choice - 1;
        print_string_matches(output, names[selected], &results[selected]);
    }
    for (i = 0; i < 4; ++i) {
        string_match_result_destroy(&results[i]);
    }
    fasta_sequences_destroy(sequences, count);
    return 1;

failure:
    fputs("String matching failed.\n", output);
    for (i = 0; i < 4; ++i) {
        string_match_result_destroy(&results[i]);
    }
    fasta_sequences_destroy(sequences, count);
    return 1;
}

static void print_benchmark_border(FILE *output,
                                   const size_t *column_widths,
                                   size_t column_count)
{
    size_t column;

    fputc('+', output);
    for (column = 0; column < column_count; ++column) {
        size_t character;

        for (character = 0; character < column_widths[column] + 2;
             ++character) {
            fputc('-', output);
        }
        fputc('+', output);
    }
    fputc('\n', output);
}

static void print_benchmark_records(FILE *output,
                                    const BenchmarkRecord *records,
                                    size_t count)
{
    static const size_t performance_widths[] = {28, 5, 13, 3, 5};
    static const size_t result_widths[] = {28, 11, 24};
    static const size_t measurement_widths[] = {28, 5, 10, 10, 10};
    size_t i;

    fputs("BENCHMARK CONFIGURATION\n", output);
    fputs("Measured records : ", output);
    fprintf(output, "%zu\n", count);
    fputs("Time unit        : milliseconds (average elapsed time)\n\n",
          output);

    fputs("ALGORITHM PERFORMANCE\n", output);
    for (i = 0; i < count; ++i) {
        size_t j;
        const char *dataset = records[i].dataset[0] == '\0'
                                  ? "-"
                                  : records[i].dataset;

        for (j = 0; j < i; ++j) {
            if (strcmp(records[j].dataset, records[i].dataset) == 0) {
                break;
            }
        }
        if (j != i) {
            continue;
        }
        fprintf(output, "\nDataset: %s\n", dataset);
        print_benchmark_border(output, performance_widths,
                               sizeof(performance_widths) /
                                   sizeof(performance_widths[0]));
        fprintf(output, "| %-28s | %5s | %13s | %3s | %-5s |\n",
                "Algorithm", "Input", "Average (ms)", "Run", "Check");
        print_benchmark_border(output, performance_widths,
                               sizeof(performance_widths) /
                                   sizeof(performance_widths[0]));
        for (j = 0; j < count; ++j) {
            if (strcmp(records[j].dataset, records[i].dataset) != 0) {
                continue;
            }
            fprintf(output,
                    "| %-28.28s | %5zu | %13.6f | %3zu | %-5s |\n",
                    records[j].algorithm, records[j].input_size,
                    records[j].average_seconds * 1000.0,
                    records[j].repetitions,
                    records[j].verified ? "PASS" : "FAIL");
        }
        print_benchmark_border(output, performance_widths,
                               sizeof(performance_widths) /
                                   sizeof(performance_widths[0]));
    }

    fputs("\nRESULTS AND COUNTS\n", output);
    print_benchmark_border(output, result_widths,
                           sizeof(result_widths) / sizeof(result_widths[0]));
    fprintf(output, "| %-28s | %-11s | %-24s |\n", "Algorithm",
            "Dataset", "Result");
    print_benchmark_border(output, result_widths,
                           sizeof(result_widths) / sizeof(result_widths[0]));
    for (i = 0; i < count; ++i) {
        char result[64];

        if (records[i].has_result_count) {
            snprintf(result, sizeof(result), "matches=%zu",
                     records[i].result_count);
        } else if (records[i].has_solution_value) {
            snprintf(result, sizeof(result), "value=%" PRId64,
                     records[i].solution_value);
        } else {
            snprintf(result, sizeof(result), "n/a");
        }
        fprintf(output, "| %-28.28s | %-11.11s | %-24.24s |\n",
                records[i].algorithm,
                records[i].dataset[0] == '\0' ? "-" : records[i].dataset,
                result);
    }
    print_benchmark_border(output, result_widths,
                           sizeof(result_widths) / sizeof(result_widths[0]));

    fputs("\nMEMORY AND OPERATION MEASUREMENTS\n", output);
    print_benchmark_border(output, measurement_widths,
                           sizeof(measurement_widths) /
                               sizeof(measurement_widths[0]));
    fprintf(output, "| %-28s | %5s | %10s | %10s | %10s |\n",
            "Algorithm", "Input", "Result B", "Aux B", "Ops");
    print_benchmark_border(output, measurement_widths,
                           sizeof(measurement_widths) /
                               sizeof(measurement_widths[0]));
    for (i = 0; i < count; ++i) {
        char result_bytes[32];
        char auxiliary_bytes[32];
        char operations[32];

        if (records[i].has_measured_result_bytes) {
            snprintf(result_bytes, sizeof(result_bytes), "%zu",
                     records[i].measured_result_bytes);
        } else {
            snprintf(result_bytes, sizeof(result_bytes), "n/a");
        }
        if (records[i].has_estimated_auxiliary_bytes) {
            snprintf(auxiliary_bytes, sizeof(auxiliary_bytes), "%zu",
                     records[i].estimated_auxiliary_bytes);
        } else {
            snprintf(auxiliary_bytes, sizeof(auxiliary_bytes), "n/a");
        }
        if (records[i].has_operation_count) {
            snprintf(operations, sizeof(operations), "%llu",
                     (unsigned long long)records[i].operation_count);
        } else {
            snprintf(operations, sizeof(operations), "n/a");
        }
        fprintf(output, "| %-28.28s | %5zu | %10s | %10s | %10s |\n",
                records[i].algorithm, records[i].input_size,
                result_bytes, auxiliary_bytes, operations);
    }
    print_benchmark_border(output, measurement_widths,
                           sizeof(measurement_widths) /
                               sizeof(measurement_widths[0]));
}

static void print_string_benchmark_summaries(
    FILE *output, const BenchmarkRecord *records, size_t count)
{
    size_t i;

    fputs("\nPERFORMANCE SUMMARY\n", output);
    for (i = 0; i < count; ++i) {
        size_t j;
        size_t algorithm_count = 0;
        size_t occurrence_count = 0;
        size_t fastest = i;
        size_t slowest = i;
        double total_seconds = 0.0;
        int occurrences_match = 1;
        int all_verified = 1;

        if (!records[i].has_result_count) {
            continue;
        }
        for (j = 0; j < i; ++j) {
            if (records[j].has_result_count &&
                records[j].input_size == records[i].input_size &&
                strcmp(records[j].dataset, records[i].dataset) == 0) {
                break;
            }
        }
        if (j != i) {
            continue;
        }
        for (j = i; j < count; ++j) {
            if (!records[j].has_result_count ||
                records[j].input_size != records[i].input_size ||
                strcmp(records[j].dataset, records[i].dataset) != 0) {
                continue;
            }
            if (algorithm_count == 0) {
                occurrence_count = records[j].result_count;
                fastest = j;
                slowest = j;
            } else {
                if (records[j].result_count != occurrence_count) {
                    occurrences_match = 0;
                }
                if (records[j].average_seconds <
                    records[fastest].average_seconds) {
                    fastest = j;
                }
                if (records[j].average_seconds >
                    records[slowest].average_seconds) {
                    slowest = j;
                }
            }
            if (!records[j].verified) {
                all_verified = 0;
            }
            total_seconds += records[j].average_seconds;
            ++algorithm_count;
        }
        if (algorithm_count == 0) {
            continue;
        }

        fprintf(output, "\nDataset: %s (input size %zu)\n",
                records[i].dataset, records[i].input_size);
        fprintf(output, "  Fastest Algorithm : %s (%.6f ms)\n",
                records[fastest].algorithm,
                records[fastest].average_seconds * 1000.0);
        fprintf(output, "  Slowest Algorithm : %s (%.6f ms)\n",
                records[slowest].algorithm,
                records[slowest].average_seconds * 1000.0);
        fprintf(output, "  Average Time      : %.6f ms\n",
                total_seconds * 1000.0 / (double)algorithm_count);
        fprintf(output,
                "  Occurrences       : %zu across %zu algorithm(s)\n",
                occurrence_count, algorithm_count);
        fprintf(output, "  Correctness Check : %s\n",
                all_verified && occurrences_match ? "PASSED" : "FAILED");
    }
}

static int run_benchmark_menu(FILE *output)
{
    BenchmarkRecord records[BENCHMARK_SUITE_MAX_RECORDS];
    size_t count = 0;

    if (!benchmark_run_suite(records, BENCHMARK_SUITE_MAX_RECORDS, &count)) {
        fprintf(output, "Benchmark suite failed after %zu record(s).\n",
                count);
        return 1;
    }
    print_benchmark_records(output, records, count);
    print_string_benchmark_summaries(output, records, count);
    fputs("Times are measurements for this run, not complexity guarantees.\n",
          output);
    return 1;
}

int genome_cli_run_benchmarks(FILE *output)
{
    if (output == NULL) {
        return 0;
    }
    return run_benchmark_menu(output);
}

typedef int (*CliWorkflow)(FILE *input, FILE *output);

static int run_module_workflow(FILE *input, FILE *output, size_t module)
{
    switch (module) {
    case 1:
        return run_dna_statistics(input, output);
    case 2:
        return run_fasta_display(input, output);
    case 3:
        return run_sort_and_search(input, output);
    case 4:
        return run_divide_and_conquer(output);
    case 5:
        return run_greedy(output);
    case 6:
        return run_lcs_workflow(input, output);
    case 7:
        return run_branch_and_bound(input, output);
    case 8:
        return run_graph_flow(output);
    case 9:
        return run_string_search(input, output);
    case 10:
        return run_benchmark_menu(output);
    default:
        return 0;
    }
}

static int run_simple_submenu(FILE *input,
                              FILE *output,
                              const CliContext *context,
                              size_t module,
                              const char *title,
                              const char *description,
                              const char *action_label)
{
    int first_display = 1;

    for (;;) {
        size_t choice;
        int status;

        if (first_display) {
            clear_screen(context, output);
            first_display = 0;
        }
        print_heading(output, context, title);
        fprintf(output, "%s\n\n[01] %s\n[00] Return to Main Menu\n\n",
                description, action_label);
        status = read_menu_choice(input, output, 1, &choice);
        if (status == 0) {
            return 1;
        }
        if (status < 0) {
            print_error(output, context,
                        "Invalid submenu choice. Choose 01 to run the module "
                        "or 00 to return.");
            continue;
        }
        if (choice == 0) {
            return 1;
        }
        (void)run_module_workflow(input, output, module);
    }
}

static void print_main_menu(FILE *output, const CliContext *context)
{
    print_heading(output, context,
                  "GENOME SEQUENCE SEARCH AND ALIGNMENT TOOLKIT");
    if (context->color_enabled) {
        fputs("\033[1;33m", output);
    }
    fputs("PROJECT MODULES\n", output);
    if (context->color_enabled) {
        fputs("\033[0m", output);
    }
    fputs("[01] DNA Sequence Operations\n"
          "[02] FASTA File Operations\n"
          "[03] Sorting and Searching Algorithms\n"
          "[04] Divide and Conquer Algorithms\n"
          "[05] Greedy Algorithms\n"
          "[06] Dynamic Programming Algorithms\n"
          "[07] Branch and Bound Algorithms\n"
          "[08] Graph and Flow Algorithms\n"
          "[09] String Matching Algorithms\n"
          "[10] Benchmarking and Performance Analysis\n"
          "[11] DNA Sequence Alignment\n\n"
          "[00] Exit Application\n\n",
          output);
}

static int confirm_exit(FILE *input,
                        FILE *output,
                        const CliContext *context)
{
    char answer[32];

    for (;;) {
        int status = read_line(input, output,
                               "Exit the application? [y/N]: ",
                               answer, sizeof(answer));
        if (status == 0) {
            return 1;
        }
        if (status < 0) {
            continue;
        }
        if (answer[0] == '\0' || strcmp(answer, "n") == 0 ||
            strcmp(answer, "N") == 0 || strcmp(answer, "no") == 0 ||
            strcmp(answer, "NO") == 0) {
            fputs("Exit cancelled.\n", output);
            return 0;
        }
        if (strcmp(answer, "y") == 0 || strcmp(answer, "Y") == 0 ||
            strcmp(answer, "yes") == 0 || strcmp(answer, "YES") == 0) {
            print_success(output, context, "Goodbye.");
            return 1;
        }
        print_error(output, context,
                    "Enter Y to exit or N to return to the main menu.");
    }
}

int genome_cli_run(FILE *input, FILE *output)
{
    CliContext context = {0};

    if (input == NULL || output == NULL) {
        return 0;
    }
    context.match_score = 2;
    context.mismatch_penalty = -1;
    context.gap_penalty = -2;
    configure_terminal(&context, output);
    for (;;) {
        size_t choice;
        int status;

        clear_screen(&context, output);
        print_main_menu(output, &context);
        status = read_menu_choice(input, output, 11, &choice);
        if (status == 0) {
            fputs("Input ended; exiting safely.\n", output);
            restore_terminal(&context);
            return 1;
        }
        if (status < 0) {
            print_error(output, &context,
                        "Invalid menu choice. Choose a module from 01 "
                        "through 11, or 00 to exit.");
            continue;
        }
        if (choice == 0) {
            if (confirm_exit(input, output, &context)) {
                restore_terminal(&context);
                return 1;
            }
            continue;
        }
        if (choice == 11) {
            (void)run_dna_alignment(input, output, &context);
            continue;
        }
        switch (choice) {
        case 1:
            (void)run_simple_submenu(
                input, output, &context, choice, "DNA SEQUENCE OPERATIONS",
                "Inspect DNA sequence composition and nucleotide statistics.",
                "Analyze a DNA sequence");
            break;
        case 2:
            (void)run_simple_submenu(
                input, output, &context, choice, "FASTA FILE OPERATIONS",
                "Load FASTA records and inspect their sequence statistics.",
                "Load and display a FASTA file");
            break;
        case 3:
            (void)run_simple_submenu(
                input, output, &context, choice,
                "SORTING AND SEARCHING ALGORITHMS",
                "Sort FASTA records and search by the selected key.",
                "Sort records and search for a sequence");
            break;
        case 4:
            (void)run_simple_submenu(
                input, output, &context, choice,
                "DIVIDE AND CONQUER ALGORITHMS",
                "Explore divide-and-conquer algorithm demonstrations.",
                "Run the maximum-subarray demonstration");
            break;
        case 5:
            (void)run_simple_submenu(
                input, output, &context, choice, "GREEDY ALGORITHMS",
                "Explore the toolkit's greedy algorithm demonstration.",
                "Run activity selection");
            break;
        case 6:
            (void)run_simple_submenu(
                input, output, &context, choice,
                "DYNAMIC PROGRAMMING ALGORITHMS",
                "Compare DNA strings using a common subsequence.",
                "Run the LCS comparison");
            break;
        case 7:
            (void)run_simple_submenu(
                input, output, &context, choice,
                "BRANCH AND BOUND ALGORITHMS",
                "Compare exact 0/1 Knapsack search strategies.",
                "Compare Knapsack solvers");
            break;
        case 8:
            (void)run_simple_submenu(
                input, output, &context, choice, "GRAPH AND FLOW ALGORITHMS",
                "Explore maximum-flow and bipartite-matching examples.",
                "Run graph and flow demonstrations");
            break;
        case 9:
            (void)run_simple_submenu(
                input, output, &context, choice,
                "STRING MATCHING ALGORITHMS",
                "Search DNA patterns with the supported string matchers.",
                "Search a FASTA sequence for a pattern");
            break;
        case 10:
            (void)run_simple_submenu(
                input, output, &context, choice,
                "BENCHMARKING AND PERFORMANCE ANALYSIS",
                "Run verified benchmarks and inspect measured performance.",
                "Run the benchmark suite");
            break;
        default:
            break;
        }
    }
}

static int run_normal_demo(FILE *output)
{
    DnaSequence *sequences = NULL;
    size_t count = 0;
    StringMatchResult results[4];
    const char *names[] = {
        "Naive", "Rabin-Karp", "KMP", "Finite Automaton"
    };
    const char *pattern = "ACG";
    size_t i;
    int successful = 0;

    if (!fasta_read_file("data/samples/demo.fasta", &sequences, &count)) {
        fputs("Normal demo could not open data/samples/demo.fasta.\n", output);
        return 0;
    }
    fprintf(output, "Scenario 1 - Normal Operation: loaded %zu FASTA records.\n",
            count);
    for (i = 0; i < count; ++i) {
        print_statistics(output, sequences[i].identifier,
                         sequences[i].bases, sequences[i].length);
    }
    if (count == 0) {
        goto cleanup;
    }

    if (!sequence_merge_sort(sequences, count, SEQUENCE_SORT_BY_LENGTH,
                             NULL)) {
        goto cleanup;
    }
    fputs("Sorted by sequence length:\n", output);
    for (i = 0; i < count; ++i) {
        fprintf(output, "  %s: %zu bases\n",
                sequences[i].identifier, sequences[i].length);
    }

    for (i = 0; i < 4; ++i) {
        string_match_result_init(&results[i]);
    }
    {
        size_t source_index;
        for (source_index = 0; source_index < count; ++source_index) {
            if (strcmp(sequences[source_index].identifier,
                       "sample_alpha") == 0) {
                break;
            }
        }
        if (source_index == count ||
            !string_match_naive(sequences[source_index].bases,
                                sequences[source_index].length,
                                pattern, strlen(pattern), &results[0]) ||
            !string_match_rabin_karp(sequences[source_index].bases,
                                     sequences[source_index].length,
                                     pattern, strlen(pattern), &results[1]) ||
            !string_match_kmp(sequences[source_index].bases,
                              sequences[source_index].length,
                              pattern, strlen(pattern), &results[2]) ||
            !string_match_finite_automaton(
                sequences[source_index].bases,
                sequences[source_index].length, pattern, strlen(pattern),
                &results[3])) {
            for (i = 0; i < 4; ++i) {
                string_match_result_destroy(&results[i]);
            }
            goto cleanup;
        }
    }
    for (i = 0; i < 4; ++i) {
        if (!results_are_equal(&results[0], &results[i])) {
            size_t result_index;
            for (result_index = 0; result_index < 4; ++result_index) {
                string_match_result_destroy(&results[result_index]);
            }
            goto cleanup;
        }
    }
    for (i = 0; i < 4; ++i) {
        print_string_matches(output, names[i], &results[i]);
        string_match_result_destroy(&results[i]);
    }
    successful = 1;
    fputs("All four string matchers agree.\n", output);

cleanup:
    fasta_sequences_destroy(sequences, count);
    if (!successful) {
        fputs("Normal demo could not complete one or more operations.\n",
              output);
    }
    return successful;
}

static int run_constraint_shock_demo(FILE *output)
{
    const size_t weights[] = {4, 6, 8, 10, 12, 14};
    const int64_t values[] = {7, 10, 13, 16, 19, 21};
    const size_t capacities[] = {30, 12};
    size_t i;

    fputs("Scenario 2 - Capacity Shock\n"
          "Items (weight:value): 4:7 6:10 8:13 10:16 12:19 14:21\n",
          output);
    for (i = 0; i < 2; ++i) {
        int64_t maximum_value = 0;
        if (!dp_knapsack_01(weights, values, 6, capacities[i],
                            &maximum_value)) {
            fputs("Knapsack demonstration failed.\n", output);
            return 0;
        }
        fprintf(output, "Capacity %zu -> maximum value %" PRId64 "\n",
                capacities[i], maximum_value);
    }
    fputs("The reduced capacity admits fewer combinations; these values were "
          "computed by the DP solver.\n", output);
    return 1;
}

typedef struct {
    const char *text;
    size_t text_length;
    const char *pattern;
    size_t pattern_length;
    const StringMatchResult *reference;
    StringMatchResult result;
    int algorithm;
} StressSearchContext;

static int stress_search_prepare(void *context_pointer)
{
    StressSearchContext *context =
        (StressSearchContext *)context_pointer;
    string_match_result_destroy(&context->result);
    return 1;
}

static int stress_search_operation(void *context_pointer)
{
    StressSearchContext *context =
        (StressSearchContext *)context_pointer;

    if (context->algorithm == 0) {
        return string_match_naive(context->text, context->text_length,
                                  context->pattern, context->pattern_length,
                                  &context->result);
    }
    if (context->algorithm == 1) {
        return string_match_rabin_karp(
            context->text, context->text_length,
            context->pattern, context->pattern_length, &context->result);
    }
    if (context->algorithm == 2) {
        return string_match_kmp(context->text, context->text_length,
                                context->pattern, context->pattern_length,
                                &context->result);
    }
    return string_match_finite_automaton(
        context->text, context->text_length, context->pattern,
        context->pattern_length, &context->result);
}

static int stress_search_verify(void *context_pointer)
{
    StressSearchContext *context =
        (StressSearchContext *)context_pointer;
    return results_are_equal(&context->result, context->reference);
}

static void stress_search_cleanup(void *context_pointer)
{
    StressSearchContext *context =
        (StressSearchContext *)context_pointer;
    string_match_result_destroy(&context->result);
}

static int run_stress_demo(FILE *output)
{
    static const char *algorithm_names[] = {
        "Naive String Matching", "Rabin-Karp", "Knuth-Morris-Pratt",
        "Finite Automaton"
    };
    static const size_t pattern_lengths[] = {12, 16, 20};
    char *text = NULL;
    DnaSequence *sequences = NULL;
    char (*identifiers)[24] = NULL;
    char (*bases)[STRESS_SEQUENCE_LENGTH + 1] = NULL;
    uint32_t random_state = UINT32_C(0xC17DAA);
    size_t i;
    size_t j;
    struct timespec start_time;
    struct timespec end_time;
    int successful = 0;

    text = (char *)malloc(STRESS_TEXT_LENGTH);
    sequences = (DnaSequence *)calloc(STRESS_SEQUENCE_COUNT,
                                      sizeof(*sequences));
    identifiers = (char (*)[24])calloc(STRESS_SEQUENCE_COUNT,
                                        sizeof(*identifiers));
    bases = (char (*)[STRESS_SEQUENCE_LENGTH + 1])malloc(
        STRESS_SEQUENCE_COUNT * sizeof(*bases));
    if (text == NULL || sequences == NULL || identifiers == NULL ||
        bases == NULL) {
        goto cleanup;
    }
    for (i = 0; i < STRESS_TEXT_LENGTH; ++i) {
        random_state = random_state * UINT32_C(1664525) +
                       UINT32_C(1013904223);
        text[i] = "ACGT"[(random_state >> 30) & UINT32_C(3)];
    }
    fputs("Scenario 3 - Large-Data Stress Test\n", output);
    fprintf(output,
            "Generated %d DNA bytes and %d sequences of %d bases each "
            "(fixed LCG seed 0x00C17DAA).\n",
            STRESS_TEXT_LENGTH, STRESS_SEQUENCE_COUNT,
            STRESS_SEQUENCE_LENGTH);
    fprintf(output,
            "Estimated input and sequence payload (excluding matcher "
            "workspace and allocator overhead): %zu bytes.\n",
            (size_t)STRESS_TEXT_LENGTH +
                STRESS_SEQUENCE_COUNT * sizeof(*sequences) +
                STRESS_SEQUENCE_COUNT * sizeof(*identifiers) +
                STRESS_SEQUENCE_COUNT * sizeof(*bases));

    for (i = 0; i < 3; ++i) {
        size_t pattern_length = pattern_lengths[i];
        size_t start = (i + 1) * (STRESS_TEXT_LENGTH / 4);
        char pattern[21];
        StringMatchResult reference;
        size_t algorithm_index;

        memcpy(pattern, text + start, pattern_length);
        pattern[pattern_length] = '\0';
        string_match_result_init(&reference);
        if (!string_match_naive(text, STRESS_TEXT_LENGTH, pattern,
                                pattern_length, &reference)) {
            string_match_result_destroy(&reference);
            goto cleanup;
        }
        fprintf(output, "Pattern length %zu: %zu occurrence(s).\n",
                pattern_length, reference.count);

        for (algorithm_index = 0; algorithm_index < 4; ++algorithm_index) {
            StressSearchContext context;
            BenchmarkRecord record;
            char record_name[BENCHMARK_ALGORITHM_NAME_CAPACITY];
            int written;

            memset(&context, 0, sizeof(context));
            context.text = text;
            context.text_length = STRESS_TEXT_LENGTH;
            context.pattern = pattern;
            context.pattern_length = pattern_length;
            context.reference = &reference;
            context.algorithm = (int)algorithm_index;
            string_match_result_init(&context.result);
            written = snprintf(record_name, sizeof(record_name),
                               "%s (%zu-byte pattern)",
                               algorithm_names[algorithm_index],
                               pattern_length);
            if (written < 0 ||
                (size_t)written >= sizeof(record_name) ||
                !benchmark_measure(record_name, STRESS_TEXT_LENGTH, 1,
                                   stress_search_prepare,
                                   stress_search_operation,
                                   stress_search_verify,
                                   stress_search_cleanup, NULL, &context,
                                   &record)) {
                stress_search_cleanup(&context);
                string_match_result_destroy(&reference);
                goto cleanup;
            }
            fprintf(output, "%s: %.6f seconds, verified.\n",
                    record.algorithm, record.average_seconds);
        }
        string_match_result_destroy(&reference);
    }

    for (i = 0; i < STRESS_SEQUENCE_COUNT; ++i) {
        (void)snprintf(identifiers[i], sizeof(identifiers[i]),
                       "synthetic_%04lu", (unsigned long)i);
        for (j = 0; j < STRESS_SEQUENCE_LENGTH; ++j) {
            random_state = random_state * UINT32_C(1664525) +
                           UINT32_C(1013904223);
            bases[i][j] = "ACGT"[(random_state >> 30) & UINT32_C(3)];
        }
        bases[i][STRESS_SEQUENCE_LENGTH] = '\0';
        sequences[i].identifier = identifiers[i];
        sequences[i].bases = bases[i];
        sequences[i].length = STRESS_SEQUENCE_LENGTH;
    }
    if (timespec_get(&start_time, TIME_UTC) != TIME_UTC ||
        !sequence_quick_sort(sequences, STRESS_SEQUENCE_COUNT,
                             SEQUENCE_SORT_BY_LENGTH, NULL) ||
        timespec_get(&end_time, TIME_UTC) != TIME_UTC) {
        goto cleanup;
    }
    for (i = 1; i < STRESS_SEQUENCE_COUNT; ++i) {
        if (sequences[i - 1].length > sequences[i].length) {
            goto cleanup;
        }
    }
    fprintf(output,
            "Quick-sorted %d records by length in approximately %.6f "
            "seconds.\n",
            STRESS_SEQUENCE_COUNT,
            difftime(end_time.tv_sec, start_time.tv_sec) +
                ((double)end_time.tv_nsec - (double)start_time.tv_nsec) /
                    1000000000.0);
    successful = 1;
    fputs("All twelve string-search benchmark records were verified.\n",
          output);

cleanup:
    free(bases);
    free(identifiers);
    free(sequences);
    free(text);
    if (!successful) {
        fputs("Stress demo failed; reduce the workload or check available "
              "memory.\n",
              output);
    }
    return successful;
}

int genome_cli_run_demo(FILE *output, const char *scenario)
{
    if (output == NULL || scenario == NULL) {
        return 0;
    }
    if (strcmp(scenario, "normal") == 0) {
        return run_normal_demo(output);
    }
    if (strcmp(scenario, "shock") == 0) {
        return run_constraint_shock_demo(output);
    }
    if (strcmp(scenario, "stress") == 0) {
        return run_stress_demo(output);
    }
    fputs("Unknown demo scenario. Choose normal, shock, or stress.\n",
          output);
    return 0;
}
