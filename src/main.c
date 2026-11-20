#include "genome/cli/cli.h"
#include "genome/io/fasta_reader.h"

#include <stdio.h>
#include <string.h>

static int print_fasta_file(const char *path)
{
    DnaSequence *sequences = NULL;
    size_t count = 0;
    size_t i;

    if (!fasta_read_file(path, &sequences, &count)) {
        fprintf(stderr, "Error: could not read valid FASTA data from '%s'.\n",
                path);
        return 1;
    }
    printf("Loaded %zu DNA sequence(s).\n", count);
    for (i = 0; i < count; ++i) {
        printf("%zu. %s | length: %zu | sequence: %s\n",
               i + 1, sequences[i].identifier, sequences[i].length,
               sequences[i].bases);
    }
    fasta_sequences_destroy(sequences, count);
    return 0;
}

int main(int argc, char **argv)
{
    if (argc == 1) {
        return genome_cli_run(stdin, stdout) ? 0 : 1;
    }
    if (argc == 2 && strcmp(argv[1], "--help") == 0) {
        puts("Usage:\n"
             "  genome_tool                 Start interactive menu\n"
             "  genome_tool <input.fasta>   Load and display FASTA records\n"
             "  genome_tool --demo <name>  Run normal, shock, or stress demo\n"
             "  genome_tool --benchmark   Run the verified benchmark suite");
        return 0;
    }
    if (argc == 3 && strcmp(argv[1], "--demo") == 0) {
        return genome_cli_run_demo(stdout, argv[2]) ? 0 : 1;
    }
    if (argc == 2 && strcmp(argv[1], "--benchmark") == 0) {
        return genome_cli_run_benchmarks(stdout) ? 0 : 1;
    }
    if (argc == 2 && argv[1][0] != '-') {
        return print_fasta_file(argv[1]);
    }
    fprintf(stderr,
            "Invalid arguments. Use --help for usage information.\n");
    return 1;
}
