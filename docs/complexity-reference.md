# Algorithm Complexity Reference

This reference summarizes the algorithms implemented by the toolkit. `n`
denotes the relevant input length, `V` vertices, `E` edges, `k` symbols or
congruences, `W` knapsack capacity, and `C` a text-pattern comparison cost.
Output storage is excluded from auxiliary-space columns unless noted. Exact
input representations and important limitations are described in
[`architecture.md`](architecture.md).

## Sequence and divide-and-conquer algorithms

| Algorithm | Time | Auxiliary space | Notes |
|---|---:|---:|---|
| FASTA parsing | O(B) | O(B) | `B` is total input bytes; storage includes returned sequences. |
| Sequence statistics | O(n) | O(1) | One pass; ambiguous bytes are counted. |
| Merge sort | O(n log n) | O(n) | Stable, bottom-up. |
| Three-way quick sort | O(n log n) average; O(n²) worst | O(log n) stack | In-place record sorting; recurses on the smaller partition. |
| Iterative binary search | O(log n) worst | O(1) | Input must be sorted by the selected key. |
| Recursive binary search | O(log n) worst | O(log n) stack | Input must be sorted by the selected key. |
| Maximum subarray | O(n log n) | O(log n) stack | Divide-and-conquer recurrence. |
| Exponentiation by squaring | O(log \|e\|) | O(1) | `e` is the exponent. |
| Celebrity identification | O(n) | O(1) | Reads an already-available relationship matrix. |
| Classical matrix multiplication | O(rkp) | O(1) workspace | Produces an `r × p` matrix; output storage is O(rp). |
| Strassen matrix multiplication | O(n^log₂7) ≈ O(n^2.807) | O(n²) | Square inputs; padded dimensions remain within a factor of two. |

## Greedy and graph algorithms

| Algorithm | Time | Auxiliary space | Notes |
|---|---:|---:|---|
| Activity selection | O(n log n) | O(n) | Sorts a copy of the input. |
| Fractional knapsack | O(n log n) | O(n) | Sorts items by value density. |
| Huffman tree construction | O(k log k) | O(k) | Min-heap and tree nodes; encoding/decoding additionally scale with input/output length. |
| Dijkstra (adjacency matrix) | O(V²) | O(V) | Requires finite nonnegative edge weights. |
| Prim (adjacency matrix) | O(V²) | O(V) | Returns a minimum spanning forest. |
| Kruskal | O(E log E) | O(V + E) | Edge sorting dominates disjoint-set processing. |
| Ford–Fulkerson with Edmonds–Karp BFS | O(VE²) | O(V + E) | Adjacency-list residual graph; parallel arcs are supported. |
| Maximum bipartite matching | O(VE) | O(V + E) | DFS augmenting paths; `V` includes both partitions. |

## Dynamic programming and search

| Algorithm | Time | Auxiliary space | Notes |
|---|---:|---:|---|
| 0/1 Knapsack, value only | O(n(W + 1)) | O(W + 1) | Pseudopolynomial in numeric capacity. |
| 0/1 Knapsack with reconstruction | O(n(W + 1)) | O(W + 1) values + O(n(W + 1)) decision bytes | Tie decisions deterministically skip the current item. |
| Longest common subsequence | O(mn) | O(mn) | Returns one LCS. |
| Needleman–Wunsch global alignment | O(mn) | O(mn) | Aligns both complete sequences; returns alignment strings in O(m+n) output space. |
| Smith–Waterman local alignment | O(mn) | O(mn) | Returns a highest-scoring local region; output space is O(m+n). |
| Matrix-chain multiplication | O(n³) | O(n²) | `n` matrices. |
| Optimal BST (successful searches) | O(n³) | O(n²) | The implementation examines all root candidates. |
| Bellman–Ford | O(VE) | O(V) | Single-source; supports signed edge weights. |
| Floyd–Warshall | O(V³) | O(V²) | All-pairs paths and negative-cycle impact. |
| Directed bitmask TSP | O(n²2ⁿ) | O(n2ⁿ) | Exact; default limit is 18 vertices, configurable up to 20. |
| Sum of Subsets | O(2ⁿ) worst case | O(n) search stack, excluding output | May return exponentially many subsets. |
| Backtracking 0/1 Knapsack | O(2ⁿ) worst case | O(n) | Exact exhaustive search. |
| FIFO/LIFO/least-cost Branch and Bound Knapsack | O(n log n + 2ⁿ) worst case | O(2ⁿ) | Exact; pruning improves many but not all inputs. |

## String matching and output

| Algorithm | Time | Auxiliary space | Notes |
|---|---:|---:|---|
| Naive string matching | O(nm) worst case | O(1) | Excludes returned positions. |
| Rabin–Karp | O(n + m) expected; O(nm) worst case | O(1) | Verifies bytes after hash equality; excludes returned positions. |
| Knuth–Morris–Pratt | O(n + m) | O(m) | LPS table; excludes returned positions. |
| Finite-automaton string matching | O(256m + n) | O(256m) | Byte alphabet transition table; excludes returned positions. |
| Benchmark harness | Depends on measured operation | Depends on test data and result | Preparation, verification, and cleanup are outside the timed interval. |

All reported bounds describe algorithm work, not wall-clock guarantees. The
benchmark module verifies results and reports measured elapsed time for
reproducible, deterministic inputs; timings vary by machine, compiler, timer
resolution, and system load. Operation counters are only available for the
sequence sort/search benchmark callbacks where the API explicitly records
them.

## P, NP, NP-hardness, and NP-completeness

- **P** contains decision problems solvable in polynomial time with respect
  to the encoded input length.
- **NP** contains decision problems whose proposed yes-certificates can be
  checked in polynomial time. This does not mean “non-polynomial.”
- A problem is **NP-hard** if every problem in NP can be transformed to it by
  a polynomial-time many-one reduction. The problem need not itself be in NP.
- A problem is **NP-complete** when it is both in NP and NP-hard.
- Whether **P = NP** is an open question; this project does not settle it.

The decision version of 0/1 Knapsack asks whether a subset can meet a capacity
and a target value. In its standard binary-encoded form it is weakly
NP-complete; the corresponding optimization problem is weakly NP-hard. The
toolkit's O(nW) dynamic program is pseudopolynomial: it is polynomial in the
numeric capacity `W`, which may be exponential in the number of bits used to
encode `W`. The exact backtracking and Branch and Bound variants can require
exponential time in the worst case.

The project's existing requirements and architecture documents do not record
a particular polynomial reduction used by the coursework. This reference
therefore explains the definitions and the Knapsack classification, but does
not claim that a specific reduction was implemented or demonstrated.
