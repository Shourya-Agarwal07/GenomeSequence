# Requirements and scope

## Project constraints

- Implement the project in standard C17, with C extensions disabled.
- Keep the existing `genome_core`, `genome_tool`, and `genome_tests` targets.
- Use the standard C library only for this milestone.
- Compile with `-Wall -Wextra -Wpedantic` on GCC-compatible compilers.
- Preserve the existing source and include directory structure.

## DNA sequence foundation

- Accept non-empty DNA strings containing only A, C, G, and T, in either case.
- Reject null or empty sequence strings and all other characters.
- Create dynamically allocated copies of non-empty identifiers and valid DNA.
- Record the number of DNA bases, excluding the terminating null character.
- On destruction, release owned strings and reset the structure to null/zero fields.
- A sequence passed to `dna_sequence_create` must be zero-initialized or destroyed
  first; successful creation transfers ownership of its allocated strings to it.

## FASTA reader

- Read one or more records with `>` headers and multiline DNA sequence data.
- Dynamically grow line, sequence, and record storage; do not silently truncate
  long input.
- Accept LF and Windows CRLF line endings, and accept uppercase or lowercase DNA.
- Ignore blank lines. Reject sequence data before the first header, an empty
  identifier, a record without bases, invalid DNA, and an empty file.
- Return failure for missing/unreadable files and allocation/read errors.
- Initialize provided output parameters to `NULL` and zero before reading. On
  success, the caller owns the returned records and releases them with
  `fasta_sequences_destroy`.
- Release all partial allocations on failure and guard capacity arithmetic
  against `size_t` overflow.

## Implemented algorithm milestones

- Milestone 2A provides merge sort, three-way quick sort, and iterative and
  recursive binary search for DNA sequence records.
- Milestone 2B provides maximum subarray, exponentiation by squaring, celebrity
  identification, classical matrix multiplication, and Strassen multiplication.
- Milestone 3 provides activity selection, fractional knapsack, Huffman coding,
  Dijkstra shortest paths, Prim minimum spanning trees, and Kruskal minimum
  spanning trees.
- Milestone 4 provides tabulated dynamic programming solutions for 0/1
  knapsack, longest common subsequence, matrix chain multiplication with an
  optimal parenthesization, and an Optimal Binary Search Tree using
  successful-search frequencies only. Bellman–Ford shortest paths use
  dynamic-programming relaxation and support signed edge weights. Floyd–Warshall
  computes all-pairs paths in a directed adjacency matrix, and bitmask
  dynamic programming solves bounded-size directed TSP instances.
- Milestone 5 provides Sum of Subsets by backtracking, an
  independent backtracking solution to 0/1 Knapsack, and FIFO, LIFO, and
  least-cost Branch and Bound solutions to the same 0/1 Knapsack problem.
- Milestone 6 provides Edmonds–Karp maximum flow and DFS augmenting-path
  maximum bipartite matching.
- Milestone 7 provides byte-oriented Naive, Rabin–Karp, KMP, and
  finite-automaton string matching.

### DNA sequence alignment

- Implement Needleman–Wunsch global alignment and Smith–Waterman local
  alignment as reusable C17 APIs and expose them in the existing CLI.
- Accept two non-empty DNA strings containing only A, C, G, and T; accept
  either input case and reject invalid characters.
- Support configurable match, mismatch, and gap scores, defaulting to +2, -1,
  and -2. Generate gaps automatically during traceback.
- Return the alignment score, aligned sequences, match indicator, match,
  mismatch, and gap counts, percentage identity, and local input regions.
- Check DP table dimensions and score arithmetic for overflow; free all
  allocations and preserve existing results on failure.
- Use O(mn) time and O(mn) auxiliary storage for input lengths `m` and `n`.
- Explain that LCS finds a gap-free common subsequence, Needleman–Wunsch aligns
  complete inputs, and Smith–Waterman finds a high-scoring local region.

### Floyd–Warshall

- Accept a nonempty `V × V` directed matrix with separate row-major weight and
  edge-presence arrays. An absent edge must be distinguishable from a
  zero-weight edge.
- Support positive, zero, and negative edge weights, disconnected vertices,
  and self-loops.
- Return all-pairs `int64_t` distances, reachability, next hops, global
  negative-cycle status, and a per-pair indication of paths affected by a
  negative cycle.
- Reject path reconstruction for invalid indices, unreachable pairs, or pairs
  affected by negative cycles; bound reconstruction against cyclic next hops.
- Check allocation dimensions and every distance addition. On failure,
  preserve any previously valid result.
- Provide a matching result destructor and document ownership of reconstructed
  paths.
- Use O(V³) time and O(V²) auxiliary storage.

Floyd–Warshall must be distinguished from Bellman–Ford: it computes every
source/destination pair and detects negative cycles anywhere in the graph;
Bellman–Ford computes from one source and reports only cycles reachable from
that source. For graphs without negative cycles, repeated Bellman–Ford runs
must agree with each Floyd–Warshall row.

### Travelling Salesman Problem

- Solve directed, possibly asymmetric graphs using the exact bitmask state
  `dp[mask][last]`, visiting every vertex once and returning to the start.
- Accept any representable signed edge cost, including zero and negative
  values. Use the same explicit edge-presence representation as
  Floyd–Warshall.
- Return a minimum cost and an allocated closed tour, or a distinct no-tour
  status. Validate the start, graph, table dimensions, shifts, and arithmetic.
- Check all allocations and free partial state on failure. Preserve an
  existing result on errors; provide a matching result destructor.
- Limit the default instance to 18 vertices; allow compile-time configuration
  only from 1 through 20 to constrain exponential memory use.
- Use O(n²·2ⁿ) time and O(n·2ⁿ) auxiliary storage. Explain the exponential
  subset-state growth and document the limitation.

### Backtracking and Branch and Bound

- Implement Sum of Subsets for positive signed integer inputs and a
  nonnegative target. Return every index subset exactly once. Equal-valued
  items at different positions are distinct choices; target zero returns the
  empty subset.
- Implement an independent recursive backtracking 0/1 Knapsack that returns
  an optimal value and selected indices. It must not call the DP Knapsack.
- Implement FIFO, LIFO, and least-cost Branch and Bound Knapsack variants
  using a shared input and result type. Each node represents the next item
  decision; include branches exceeding capacity are infeasible.
- Sort the Branch and Bound items by descending value-to-weight ratio with
  original indices retained. Treat positive-value zero-weight items safely.
  The ratio-based upper bound must be optimistic, saturate rather than
  overflow, and only prune nodes that cannot improve the best feasible value.
- Resolve equal ratios and equal node priorities deterministically.
- Reject negative values, invalid pointers, allocation-size overflow, and
  unrepresentable feasible value sums. Free partial allocations and preserve
  an existing result on failure.
- Document ownership and provide result destructors. Selection indices must
  be unique, in range, ascending, and reconstruct the reported value within
  capacity.
- Use a test-only brute-force oracle for many small deterministic-random
  Knapsack instances. Compare the DP implementation, backtracking, and all
  three Branch and Bound methods.

Sum of Subsets has O(2ⁿ) worst-case search. Backtracking Knapsack and each
Branch and Bound method examine at most O(2ⁿ)
state-space nodes. The ratio upper bound is O(1) per node after O(n log n)
item sorting and suffix preprocessing. Backtracking uses O(n) auxiliary stack
and path storage, excluding returned solutions; Branch and Bound uses space
proportional to generated nodes. Sum of Subsets output storage is proportional
to the total number of returned indices.

### Ford–Fulkerson maximum flow

- Accept a directed edge-list network with signed `int64_t` capacities. Reject
  negative capacities and invalid source/sink vertices; require distinct
  source and sink.
- Use BFS augmenting paths (Edmonds–Karp), paired reverse residual arcs, and
  per-edge flow output. Preserve parallel input edges as separate arcs.
- Return the maximum-flow value, final source-side reachability, crossing
  minimum-cut edges, total cut capacity, and flow corresponding to each input
  edge.
- Check allocation dimensions, residual updates, total flow, and cut-capacity
  arithmetic. On invalid input, allocation failure, or unrepresentable
  maximum flow, preserve any previous result.
- Provide a result destructor. The result owns all returned arrays.
- Use O(VE²) time and O(V+E) auxiliary storage with adjacency lists and paired
  residual arcs.

### Maximum bipartite matching

- Accept left/right partition sizes and a list of left-to-right edges.
- Use DFS augmenting paths; return the matching size, matched pairs, and
  `SIZE_MAX`-marked left-to-right and right-to-left maps.
- Validate all edge indices and allocation-size calculations. Permit empty
  partitions and isolated vertices.
- Treat duplicate edges as redundant adjacency entries: they do not increase
  matching size or create duplicate pairs.
- Preserve a previous result on invalid input or allocation failure and
  provide a result destructor.
- Use O(VE) time for the basic augmenting-path search and O(V+E) storage.

### String matching

- Accept explicit-length text and pattern byte strings with `size_t` lengths.
  Support arbitrary bytes (including embedded NUL), zero-based positions,
  overlapping occurrences, and patterns longer than the text.
- Define an empty pattern as matching every position from zero through the
  text length. Require non-NULL input pointers even for zero-length strings.
- Provide Naive O(nm), Rabin–Karp O(n + m) expected / O(nm) worst-case, KMP
  O(n + m), and finite-automaton O(|Σ|m + n) implementations. Rabin–Karp must
  verify bytes after a hash match; KMP must use an LPS table. The automaton
  must support the full byte alphabet, including embedded NUL bytes.
- Return all positions through a shared caller-owned result structure.
  Successful calls replace the prior result; failed calls preserve it. Check
  allocation-size arithmetic and release partial results and KMP tables.

## Milestone 8 integration and presentation

- Preserve all existing C17 algorithms and expose them through reusable APIs;
  do not replace assessed algorithms with library implementations.
- Provide an interactive command-line menu, keep the existing FASTA-path
  invocation compatible, and support named normal, constraint-shock, and
  stress demonstrations.
- Provide sequence base counts, ambiguous-base counts, and a GC percentage
  only when the sequence contains at least one canonical base and no
  ambiguous bases.
- Provide a benchmark harness with setup and verification outside the timed
  region, deterministic inputs, correctness checks, multiple problem sizes,
  elapsed times, and operation counts where existing algorithm APIs expose
  them. Include small, medium, large, and adversarial string-matching
  datasets; compare all four matchers on exactly the same input in each
  dataset. Report input size and solution quality where applicable, plus
  memory usage with measured memory clearly distinguished from estimates.
  Timings are measurements, not hard-coded pass/fail thresholds.
- Include three reproducible demos: normal FASTA/statistics/matching
  operations; a 0/1 Knapsack capacity shock; and a large-data test that
  searches multiple patterns with all four string matchers, verifies
  equivalent results, sorts a large sequence collection, records timings and
  an approximate memory observation, and reports the measured results.
- Document time and auxiliary-space complexity for all existing algorithms.
  Explain P, NP, NP-hardness, NP-completeness, the open P-versus-NP question,
  and the pseudopolynomial nature of 0/1 Knapsack DP.
- Do not invent a project-specific polynomial reduction. If the project
  materials do not specify one, explicitly state that no particular reduction
  is documented or demonstrated.
- State that the bytewise LCS demonstration is not a biological alignment
  model. Include a short presentation script and a PowerShell runner for the
  three demos.
- Add integration coverage without removing or weakening prior test cases;
  run the complete C17 build and test suite.
