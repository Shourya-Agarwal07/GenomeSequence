#ifndef GENOME_CLI_CLI_H
#define GENOME_CLI_CLI_H

#include <stdio.h>

int genome_cli_run(FILE *input, FILE *output);

int genome_cli_run_demo(FILE *output, const char *scenario);

int genome_cli_run_benchmarks(FILE *output);

#endif
