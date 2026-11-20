# Test plan

Run the `genome_tests` executable through CTest. Each case prints an explicit
`PASS` or `FAIL` result. FASTA tests write local fixture files and remove them
after each read.

## DNA sequence tests

- Valid uppercase and lowercase DNA, mixed case, and a single nucleotide.
- Invalid symbols, empty strings, and null sequence input.
- Successful creation, copied content, and correct sequence length.
- Destruction resets all fields; repeated destruction is safe.
- Empty and null identifiers are rejected.
- Invalid sequence creation fails without leaving allocated output fields.

## FASTA reader tests

- A single record and multiple records.
- Multiline sequence concatenation.
- Invalid DNA, sequence data before the first header, and empty identifiers.
- Records without DNA and empty files.
- Missing input files and output-parameter initialization on failure.
- Windows CRLF line endings.
- DNA longer than the former fixed-size line buffer.

## Sequence algorithm tests

Each sort is exercised by length and identifier, covering empty and
single-record inputs, already sorted and reverse-sorted inputs, duplicate
values, and a deterministic shuffled order. Tests check ordering and that all
input records remain present.

Each binary search implementation is exercised with both sort criteria. Tests
cover empty and single-record arrays, matches at the first and last positions,
duplicate lengths, and missing targets. A found result must include a usable
index; a missing result sets the index to `SIZE_MAX`.

Complexity and operation-counter conventions are documented in
[`architecture.md`](architecture.md). Tests establish behavior rather than
timing or asserting hardware-dependent performance.

## Divide-and-conquer tests

- Maximum subarray checks the canonical mixed-sign case, all-negative input,
  zeroes, one element, empty/NULL input, and `int64_t` sum overflow.
- Exponentiation checks positive, zero, negative, and `INT_MIN` exponents,
  `0^0`, negative bases, NULL output, non-finite input, and overflow.
- Celebrity identification checks a known celebrity, no celebrity, one person,
  empty/NULL input, a 64-person matrix, and invalid matrix-size arithmetic.
- Classical multiplication checks a rectangular product, a zero inner
  dimension, incompatible dimensions, invalid allocation dimensions, NULL
  output construction, and floating-point overflow.
- Strassen is checked against classical multiplication for a 3 × 3 matrix and
  a 65 × 65 non-power-of-two matrix. It also checks non-square inputs and
  floating-point overflow.

The matrix tests use small integer-valued doubles for exact classical expected
values and a relative tolerance when comparing Strassen with classical
multiplication.

## Dynamic programming tests

- 0/1 knapsack checks the standard example, zero-weight items, empty input,
  zero capacity, no fitting items, deterministic tie handling, unique/in-range
  selected indices, and reconstructed weight/value totals. Invalid input and
  allocation-size overflow must fail with outputs reset.
- Longest common subsequence checks empty inputs, identical and disjoint
  strings, the standard example, multiple valid answers, single characters,
  DNA sequences, NULL inputs, a larger input pair, and table-size overflow.
  Every returned LCS must be a subsequence of both inputs, and its reported
  length must match the reconstructed string.
- DNA alignment checks Needleman–Wunsch identical and different sequences,
  unequal lengths, a gap-requiring example, and case-insensitive base
  handling. Smith–Waterman checks a known local region and the zero-score
  no-positive-region case. Invalid/empty inputs and score overflow must fail;
  failed calls preserve an existing result. A scripted menu test exercises
  the global alignment flow, default scoring, traceback strings, and metrics.
  Small expected scores are calculated directly from the selected scoring
  parameters.
- Matrix chain multiplication checks the standard three-matrix example,
  single- and two-matrix chains, deterministic tie handling, and a chain with
  an overflowing candidate but a representable optimum. Invalid dimensions,
  table-size overflow, and chains whose costs cannot be represented must fail
  with outputs reset.
- Optimal Binary Search Tree checks empty, single-key, two-key, and multi-key
  inputs; sorted-key validation; duplicate rejection; zero and large
  frequencies; prefix-sum and cost overflow; NULL arguments; and repeated
  construction/destruction. A separate exhaustive recursive reference checks
  small cases. Returned trees must contain each key once in sorted in-order
  traversal, preserve each frequency, and have a reconstructed search cost
  equal to the DP result. Allocation-size overflow is tested; direct allocator
  failure injection is unavailable in this test harness.
- Bellman–Ford checks single-vertex and disconnected graphs, known shortest
  distances, multiple shortest paths, negative and zero edges, reachable and
  unreachable negative cycles, source-to-source and ordinary path
  reconstruction, invalid destinations and edges, NULL inputs, signed
  distance overflow, allocation-size overflow, repeated successful calls,
  and failed-call result preservation. Reconstructed paths must follow
  directed edges and sum to the reported destination distance; cyclic
  predecessor chains must be rejected.
- Floyd–Warshall checks one- and two-vertex matrices, positive, zero, and
  negative edges, disconnected vertices, self paths, path reconstruction,
  reachable and disconnected negative cycles, invalid matrices, NULL inputs,
  dimension and distance overflow, repeated calls, and failure preservation.
  Zero-cost edges are tested separately from absent edges using the explicit
  presence matrix. Negative-cycle tests verify both the global flag and
  per-pair affected status.
- Bellman–Ford/Floyd–Warshall comparison runs Bellman–Ford from every vertex
  on graphs without negative cycles and compares each reachable distance and
  reachability flag against the corresponding Floyd–Warshall row.
- Bitmask TSP checks single- and two-vertex tours, known and asymmetric
  directed examples, negative costs, missing edges, no-tour outcomes, invalid
  starts, the configured maximum size, and result reuse/cleanup. Each returned
  tour is checked for its start/end vertex, uniqueness of interior vertices,
  edge existence, and total cost. Small cases are cross-checked against a
  test-only brute-force permutation reference.

## Backtracking and Branch and Bound tests

- Sum of Subsets checks empty input, target zero, no solution, one and several
  solutions, duplicate values at different indices, a single matching item,
  targets larger than the total, NULL inputs, negative targets and values,
  zero-value rejection, index uniqueness, and allocation-size overflow.
  Returned subsets are checked to sum to the target, and duplicate-valued
  items are verified as distinct index choices.
- Backtracking Knapsack and FIFO, LIFO, and least-cost Branch and Bound each
  check the standard example, empty input, zero capacity, zero-weight items,
  all items fitting, no item fitting, multiple optima, invalid pointers and
  negative values, value overflow, and failed-call result preservation.
  Returned indices must be unique, in range, ascending, within capacity, and
  sum to the reported optimum.
- A deterministic pseudo-random collection of 64 small Knapsack instances is
  checked against a test-only brute-force enumerator. Every instance compares
  the existing DP method, backtracking, and all three Branch and Bound
  methods, including reconstructed result validity.

## Flow network and matching tests

- Ford–Fulkerson/Edmonds–Karp checks the standard maximum-flow value of 23,
  no-path and zero-capacity networks, multiple paths, a direct edge, and a
  network whose maximum requires cancelling an earlier path through reverse
  residual arcs. Tests also cover parallel arcs, disconnected vertices,
  invalid endpoints, equal source/sink, negative capacities, arithmetic and
  allocation-size overflow, and result preservation after failure.
- For each returned flow, verify per-edge capacity limits, flow conservation
  at intermediate vertices, and agreement between source flow, sink flow,
  reported maximum flow, and the reported minimum-cut capacity. Every cut edge
  must cross from a source-side reachable vertex to a non-reachable vertex.
- Bipartite matching checks perfect and partial matchings, no edges, isolated
  vertices, imbalanced and empty partitions, augmenting-path reassignment,
  duplicate edges, invalid endpoints, NULL pointers, and result preservation.
  Returned pairs must be input edges, and each left and right vertex may occur
  in at most one pair.
- A deterministic set of small bipartite graphs is compared with a test-only
  brute-force matching reference.

## String-matching tests

- Naive matching checks beginning/end matches, absent and longer patterns,
  repeated characters, overlaps, empty text, empty patterns, and the defined
  empty-pattern positions from zero through text length.
- Rabin–Karp checks rolling updates across many windows and uses a deliberately
  colliding pair of byte strings to verify that equal hashes are confirmed by
  byte comparison.
- KMP checks repeated-prefix patterns including `AAAA`, `ABABAC`, and
  `ABABCABAB`, as well as overlapping and full-text matches.
- All three algorithms are cross-checked against the same expected positions
  for deterministic randomized byte inputs, including embedded NULs. A large
  input case, NULL arguments, result reuse/failure preservation, and
  length/allocation overflow guards are also checked.
- Every successful result is checked for its occurrence count, sorted
  zero-based indices, and valid capacity; all result storage is destroyed.

## Greedy algorithm tests

- Dijkstra rejects negative or non-finite matrix weights.
- Dijkstra uses nonnegative `double` weights; Bellman–Ford is the directed
  `int64_t` edge-list alternative for negative weights and reachable
  negative-cycle detection.
- Prim accepts finite negative and zero weights, produces a forest for
  disconnected undirected graphs, and rejects asymmetric or non-finite
  adjacency matrices. Prim and Kruskal total weights are compared on multiple
  connected graphs with negative weights.
- Kruskal accepts finite negative weights.
- Huffman covers normal and single-symbol round trips, empty input, unknown
  symbols, invalid bits after partial decoding, invalid frequencies, duplicate
  symbols, and rejection of alphabets larger than 256 unique byte symbols.
- Huffman trees are caller-owned and released with
  `greedy_huffman_destroy`; encode/decode buffers are caller-owned and released
  with `free`. Failed tree construction and failed decoding release all
  partial allocations and leave output pointers empty.

## Milestone 8 integration tests

- Sequence statistics verify canonical base counts, ambiguous-byte counts,
  defined and undefined GC percentage, empty input, and invalid NULL input.
- Benchmark callbacks verify that setup/verification/cleanup surround the
  timed operation, repeated operations are verified, and invalid measurement
  metadata is rejected. The complete 42-record suite is run and every record
  must be marked verified. The string-matching records cover four datasets
  with identical inputs for Naive, Rabin–Karp, KMP, and finite automaton;
  memory fields, match counts, and dataset labels are checked.
- The interactive menu is exercised with scripted input across all 11
  remaining modules: DNA
  statistics, FASTA loading, sorting/searching, divide-and-conquer, greedy,
  LCS, global and local DNA alignment, Knapsack search, flow, string matching,
  and benchmarking. The test checks representative outputs,
  custom alignment scoring, and that all four string matchers agree.
- The normal demo is checked for sequence statistics and matcher agreement;
  the capacity-shock demo is checked for both capacity cases; the stress demo
  runs the deterministic million-byte text and 512-record sort, verifies all
  twelve string-search benchmark records and prints their timings.
- The suite contains 107 test groups in total. CTest runs the test executable
  from the project root so the bundled sample FASTA path is stable.
- Wall-clock durations are displayed but are not compared to fixed
  thresholds. Benchmark memory reports distinguish measured returned-result
  buffer capacity from estimated algorithm auxiliary storage; neither is a
  process-wide peak-memory measurement. C17 has no portable RSS interface.
