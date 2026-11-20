cmake_minimum_required(VERSION 3.20)

project(GenomeToolkit VERSION 0.1.0 LANGUAGES C)

set(CMAKE_C_STANDARD 17)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(genome_core
    src/domain/dna_sequence.c
    src/io/fasta_reader.c
    src/algorithms/sequence_algorithms.c
    src/algorithms/dna_alignment.c
    src/algorithms/maximum_subarray.c
    src/algorithms/fast_exponentiation.c
    src/algorithms/celebrity.c
    src/algorithms/matrix_multiplication.c
    src/algorithms/dynamic_programming.c
    src/algorithms/backtracking_algorithms.c
    src/algorithms/greedy_algorithms.c
    src/graph/ford_fulkerson.c
    src/graph/bipartite_matching.c
    src/string_matching/string_matching.c
    src/sequence/sequence_statistics.c
    src/benchmark/benchmark.c
    src/cli/cli.c
)

target_include_directories(genome_core PUBLIC include)

if(MSVC)
    target_compile_options(genome_core PRIVATE /W4)
else()
    target_compile_options(genome_core PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(genome_tool src/main.c)
target_link_libraries(genome_tool PRIVATE genome_core)
if(MSVC)
    target_compile_options(genome_tool PRIVATE /W4)
else()
    target_compile_options(genome_tool PRIVATE -Wall -Wextra -Wpedantic)
endif()

enable_testing()
add_executable(genome_tests tests/test_main.c)
target_link_libraries(genome_tests PRIVATE genome_core)
if(MSVC)
    target_compile_options(genome_tests PRIVATE /W4)
else()
    target_compile_options(genome_tests PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME genome_tests COMMAND genome_tests)
set_tests_properties(genome_tests PROPERTIES WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
