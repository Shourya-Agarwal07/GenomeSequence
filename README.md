# Genome Sequence Search and Alignment Toolkit

A C17 implementation for the DAA Project-Based Learning assignment.

## Current milestone

This starter contains:
- A DNA sequence data structure and validation.
- A basic FASTA file reader.
- Merge sort and three-way quick sort by sequence length or identifier.
- Iterative and recursive binary search by either sorting criterion.
- Divide-and-conquer maximum subarray, fast exponentiation, and celebrity search.
- Classical rectangular and Strassen square matrix multiplication.
- Tabulated dynamic programming for 0/1 knapsack with item reconstruction and
  longest common subsequence string reconstruction, plus matrix chain
  multiplication with an optimal parenthesization and an Optimal Binary
  Search Tree using successful-search frequencies, Bellman–Ford single-source
  paths, Floyd–Warshall all-pairs paths, exact bitmask-DP Travelling Salesman
  tours, and gap-based Needleman–Wunsch global and Smith–Waterman local DNA
  alignments.
- Backtracking Sum of Subsets, plus independent backtracking,
  FIFO Branch and Bound, LIFO Branch and Bound, and least-cost Branch and Bound
  implementations of 0/1 Knapsack.
- Edmonds–Karp maximum flow with minimum-cut reporting, and maximum bipartite
  matching by DFS augmenting paths.
- Naive, Rabin–Karp, KMP, and finite-automaton matching for arbitrary byte
  strings.
- A menu-driven CLI for sequence statistics, FASTA operations, sorting/searching,
  representative algorithm demonstrations, matching, benchmarks, and
  sequence alignment. The original `genome_tool <input.fasta>` form is
  retained.
- Sequence statistics and a verified benchmark suite.
- Three deterministic demonstrations: normal operation, a knapsack capacity
  shock, and large-data stress testing.
- An automated test executable covering the domain, FASTA, all remaining
  algorithm modules, and CLI integration (107 test groups).

Complexity analysis, operation-counter definitions, algorithm pseudocode, and
design decisions are documented in [docs/architecture.md](docs/architecture.md)
and summarized in the [complexity reference](docs/complexity-reference.md).

## Floyd–Warshall and bitmask TSP

The dynamic-programming graph APIs use row-major directed adjacency matrices
with a parallel edge-presence array. This distinguishes a missing edge from a
valid zero-cost edge and supports signed `int64_t` weights. Floyd–Warshall
returns all-pairs distances and reports both global negative-cycle detection
and the pairs whose shortest paths are affected. It runs in O(V³) time and
uses O(V²) auxiliary space.

The exact directed TSP solver uses `dp[mask][last]` and returns a closed tour
from the requested start, or a distinct no-tour status. The default limit is
18 vertices; it may be configured at compile time up to 20 because the
algorithm requires O(n·2ⁿ) space and O(n²·2ⁿ) time. Graph results and
reconstructed paths/tours have explicit caller ownership; see the public
header and architecture documentation for destruction functions and details.

## Backtracking and Branch and Bound

The separate `backtracking_algorithms` module implements Sum of Subsets and
three independent exact 0/1 Knapsack search strategies. Sum of Subsets reports
original input indices; equal-valued inputs at different positions represent
distinct choices.

The Branch and Bound variants use a deterministic value-density ordering and
an overflow-safe optimistic upper bound. FIFO uses breadth-first node
selection, LIFO uses depth-first selection, and least-cost uses a priority
heap that expands the node with the highest remaining value bound. All
Knapsack search results use the same caller-owned index-array type and can be
compared directly with the Dynamic Programming implementation. Search is
exponential in the worst case; see [docs/architecture.md](docs/architecture.md)
for states, pruning, complexity, and ownership details.

## Flow networks and bipartite matching

The flow module accepts directed `int64_t`-capacity edge lists and uses BFS
augmenting paths (Edmonds–Karp), including paired reverse residual arcs so
previous choices can be rerouted. It returns per-input-edge flows and the
minimum cut as well as the maximum-flow value. Parallel arcs are preserved.
The complexity is O(VE²) time and O(V+E) auxiliary space.

The bipartite matcher accepts left-to-right edge lists and uses DFS augmenting
paths. It returns matching pairs and maps in both directions; unmatched map
entries are `SIZE_MAX`. Duplicate edges are redundant. Its basic search takes
O(VE) time and O(V+E) space.

Initialize each result to `{0}` and destroy it with its matching public
destroy function. Successful calls replace prior results; failed calls
preserve them. Result arrays are caller-owned. Build and run the complete
suite with `cmake --build build` and `ctest --test-dir build --output-on-failure`;
the [test plan](docs/test-plan.md) describes the coverage.

## String matching

The string-matching API in `include/genome/string_matching/string_matching.h`
accepts explicit-length byte strings, so embedded NUL bytes are supported.
Naive search takes O(nm) worst-case time; Rabin–Karp uses a rolling hash with
byte-by-byte collision verification and O(n + m) expected time; KMP uses an
LPS (longest proper prefix/suffix) table for O(n + m) time; finite-automaton
matching uses a byte-alphabet transition table for O(256m + n) preprocessing
and search time. All return every zero-based occurrence, including overlaps.
The empty pattern matches at every position from zero through the text length.
Initialize and destroy result structures with their matching API functions.

## Greedy algorithm input and result ownership

Fractional knapsack accepts zero-value items and rejects negative values.
Weights must be nonnegative; positive-value zero-weight items contribute their
full value without consuming capacity.

Initialize `ShortestPathResult`, `PrimResult`, and `KruskalResult` to `{0}`
before first use. A successful Dijkstra, Prim, or Kruskal call replaces and
destroys any result already held by that structure. A failed call leaves the
previous result unchanged. Destroy each result with its matching destroy
function when it is no longer needed, and do not copy owning result structures
by value.

## Requirements

- Windows 11
- Visual Studio Build Tools with the **Desktop development with C++** workload (includes the MSVC C compiler)
- CMake
- VS Code (editor)

## Build with MSVC

Open a Developer Command Prompt for Visual Studio, or initialize the environment in a VS Code terminal:

```cmd
call "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\Common7\Tools\VsDevCmd.bat" -arch=x64
```

If your Visual Studio installation uses a different version or path, use its own `VsDevCmd.bat` location.

From this project folder:

```cmd
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
```

## Run the demo

```cmd
.\build\Release\genome_tool.exe data\samples\demo.fasta
```

## Milestone 8 CLI, demos, and benchmarks

Run the interactive menu with no arguments, or start a focused, repeatable
scenario:

```powershell
.\build\Release\genome_tool.exe
.\build\Release\genome_tool.exe --demo normal
.\build\Release\genome_tool.exe --demo shock
.\build\Release\genome_tool.exe --demo stress
.\build\Release\genome_tool.exe --benchmark
.\scripts\run_demos.ps1
```

The normal scenario reads `data/samples/demo.fasta`, prints per-sequence base
counts and GC percentage, demonstrates sorting, and prints cross-validated
string-matching results. The stress scenario generates a deterministic
one-million-byte DNA text and 512 synthetic sequences, benchmarks all four
string matchers on three patterns, verifies their results against Naive
matching, sorts the synthetic sequences, and prints the verified timings.
Benchmark measurements are machine-dependent observations; see
the [benchmark notes](docs/architecture.md#benchmarking),
[test plan](docs/test-plan.md), and
[presentation script](docs/presentation-script.md).

String-matching benchmarks compare Naive, Rabin–Karp, KMP, and finite
automaton matching on identical deterministic small, medium, large, and
adversarial inputs. Each row reports correctness, elapsed time, input bytes,
and match count. It separately reports the measured result-buffer allocation
and estimated algorithm auxiliary storage; neither is a process-wide peak
memory measurement. C17 has no portable process-memory API, so no RSS or
working-set measurement is claimed. Timings and memory observations are not
fabricated or used as pass/fail thresholds.

LCS is a gap-free common-subsequence comparison: it can omit characters but
does not score substitutions or construct an alignment. Needleman–Wunsch
aligns both complete DNA strings with match, mismatch, and gap scores;
Smith–Waterman uses the same scoring model but reports the best local regions.
See the [alignment user guide](docs/user-guide.md) for CLI instructions and
examples. The byte-oriented LCS demonstration remains available in the
Dynamic Programming menu and is not a biological alignment model. The
complexity discussion does not claim that P equals or differs from NP; the
coursework documentation does not specify a particular polynomial reduction.
See the [complexity reference](docs/complexity-reference.md).

## Project structure

- `include/` public header files
- `src/domain/` core DNA data model
- `src/io/` input and output functionality
- `src/algorithms/` reusable sequence, divide-and-conquer, greedy, dynamic
  programming, backtracking, and Branch and Bound algorithms
- `src/graph/` flow-network and bipartite matching algorithms
- `src/string_matching/` byte-oriented string-search algorithms
- `src/sequence/` sequence statistics
- `src/benchmark/` timing and correctness-checked benchmark suite
- `src/cli/` interactive menu and reproducible demo scenarios
- `tests/` automated tests
- `data/` reproducible sample datasets
- `docs/` project planning and design notes
- `scripts/` environment and demo runner scripts

## Development rule

Core assessed algorithms must be implemented by the student. Standard C library facilities may be used for ordinary input/output and memory management, but must not replace the algorithms being assessed.
