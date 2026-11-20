# Architecture

- `domain`: DNA sequence model and validation.
- `io`: FASTA input handling.
- `algorithms`: sorting and searching DNA sequence records.
- `sequence`: DNA sequence statistics.
- `benchmark`: repeatable measurements with correctness verification.
- `export`: streaming CSV/JSON result export.
- `cli` and `main.c`: interactive command-line entry point and demos.
- `tests`: automated tests.

The public algorithm API is in
`include/genome/algorithms/sequence_algorithms.h`; implementations are in
`src/algorithms/sequence_algorithms.c`. Both sorting routines reorder the
`DnaSequence` structs in place (the owned identifier and bases pointers move
with each struct). Callers select either sequence length or identifier as the
key. Identifier sorting uses case-sensitive `strcmp`, matching standard C
string ordering. Merge sort is stable; quick sort is not guaranteed stable.

The sort functions accept an optional statistics structure. Comparisons count
record-key comparisons, assignments count whole-record copies made by the
sorting routine, and swaps count non-self exchanges. Character comparisons
inside `strcmp`, scalar/index assignments, and search-result output writes are
not included. Search functions require records sorted with the same criterion;
they return FOUND, NOT_FOUND, or ERROR and set the output index to `SIZE_MAX`
when no match is returned.

## Algorithm explanations and pseudocode

### DNA sequence alignment

The public alignment API is
`include/genome/algorithms/dna_alignment.h`; the implementation is
`src/algorithms/dna_alignment.c`. Both routines accept non-empty strings made
only of A, C, G, and T (case-insensitive) and caller-selected signed scores.
The CLI defaults are match `+2`, mismatch `-1`, and gap `-2`. Returned
alignment strings and the match-indicator line are owned by the result and
released by `dna_alignment_result_destroy`; initialize results to `{0}`.
Invalid inputs, allocation-size overflow, allocation failure, or score
arithmetic overflow fail without replacing an existing result.

Needleman–Wunsch globally aligns both complete inputs. Its boundary rows and
columns accumulate gap penalties, and each interior cell chooses the best of
the diagonal substitution, vertical gap, or horizontal gap:

```text
F[0,j] = j * gap
F[i,0] = i * gap
F[i,j] = max(F[i-1,j-1] + substitution,
             F[i-1,j] + gap,
             F[i,j-1] + gap)
```

Smith–Waterman uses zero-valued boundaries and also allows each interior
score to reset to zero:

```text
H[i,j] = max(0,
             H[i-1,j-1] + substitution,
             H[i-1,j] + gap,
             H[i,j-1] + gap)
```

It records the highest-scoring cell and traces back only until a zero-score
cell, reporting each matching half-open input range as a 1-based inclusive
range in the CLI. If no positive local score exists, the result is an empty
alignment with score zero. Both algorithms retain a direction per table cell
and reverse the traceback output to produce the aligned strings; gaps are
generated automatically. Match identity is matches divided by alignment
columns, including gap columns.

LCS and these alignment algorithms solve different problems. LCS maximizes the
length of an ordered common subsequence and does not model mismatches or gap
costs. Needleman–Wunsch scores a complete end-to-end alignment; Smith–Waterman
finds a high-scoring local region. For input lengths `m` and `n`, each takes
O(mn) time and O(mn) auxiliary space for scores and traceback directions.
Returned alignment storage is O(m+n).

### Merge sort

Bottom-up merge sort repeatedly merges adjacent sorted runs of doubling size.
Equal keys are taken from the left run first, preserving their original order.

```text
width = 1
while width < count:
    for each adjacent pair of runs of size width:
        merge them into temporary storage
        copy the merged run back
    double width
```

### Quick sort

Three-way partitioning around the middle record's key divides a range into
less-than, equal-to, and greater-than groups. The smaller outer partition is
sorted recursively and the larger is processed in the loop, bounding call
stack usage even for unbalanced partitions.

```text
while range has at least two records:
    copy middle key as pivot
    partition into less / equal / greater groups
    recursively sort the smaller outer group
    continue with the larger outer group
```

### Iterative binary search

The search repeatedly halves a half-open range `[low, high)` and returns the
index when the midpoint key matches.

```text
low = 0; high = count
while low < high:
    middle = low + (high - low) / 2
    compare middle key with target
    return middle if equal
    retain the half that can contain the target
return NOT_FOUND
```

### Recursive binary search

The recursive version applies the same midpoint comparison to a progressively
smaller half-open range. It reports NOT_FOUND when that range becomes empty.

```text
search(low, high):
    if low >= high: return NOT_FOUND
    middle = low + (high - low) / 2
    compare middle key with target
    return middle if equal
    recurse into the possible half
```

## Complexity

Let `n` be the number of records. Comparison counts refer to record-key
comparisons; assignment counts are whole-record assignments as defined above.

| Algorithm | Best time | Average time | Worst time | Auxiliary space |
|---|---:|---:|---:|---:|
| Merge sort | O(n log n) | O(n log n) | O(n log n) | O(n) |
| Quick sort | O(n log n) | O(n log n) | O(n²) | O(log n) call stack |
| Iterative binary search | O(1) | O(log n) | O(log n) | O(1) |
| Recursive binary search | O(1) | O(log n) | O(log n) | O(log n) call stack |

| Algorithm | Comparison count | Assignments / swaps |
|---|---|---|
| Merge sort | Θ(n log n) in best, average, and worst cases | Θ(n log n) record assignments; no swaps |
| Quick sort | Θ(n log n) best/average; O(n²) worst | O(n²) assignments and swaps worst-case; at most 1 pivot copy plus 3 record assignments per counted swap |
| Iterative binary search | 1 best; O(log n) worst | No record assignments or swaps |
| Recursive binary search | 1 best; O(log n) worst | No record assignments or swaps |

The exact search count depends on the target position and data size. Quick sort's
three-way partition handles equal keys in one partition pass, while its
worst-case time remains quadratic for sufficiently unbalanced distinct-key
partitions.

## Divide-and-conquer algorithms

The Milestone 2B public API and shared types are in
`include/genome/algorithms/divide_conquer.h`. Implementations are separated by
algorithm in `src/algorithms/`. `DcaOperationCounts` records value or matrix
relationship comparisons and arithmetic additions/multiplications; allocation,
index calculations, and copies are not counted. Counters are optional.

### Maximum subarray

The implementation recursively finds the best interval in the left half, the
right half, and the interval that crosses the midpoint. The returned ending
index is inclusive. Empty input is rejected because it has no non-empty
maximum subarray; all-negative arrays still return their least-negative
element.

```text
solve(low, high):
    if high - low == 1: return that element and its index
    split at middle
    left = solve(low, middle)
    right = solve(middle, high)
    crossing = best suffix of left half + best prefix of right half
    return the best of left, right, and crossing
```

For `n > 1`, the recurrence is `T(n) = 2T(n/2) + Θ(n)`: each level scans
each element at most a constant number of times, and there are `log₂ n` levels.
By the Master Theorem, best, average, and worst time are all Θ(n log n);
auxiliary space is Θ(log n) for the recursion stack. Comparisons and additions
are Θ(n log n). Values and sums use `int64_t`; the call fails if a sum needed
by the algorithm overflows that representation.

### Exponentiation by squaring

The exponent is reduced by half on every iteration; set bits select which
squared factors are multiplied into the result. Negative exponents use the
reciprocal of the base. `0^0` is defined as `1`; zero to a negative exponent is
rejected.

```text
result = 1; factor = base; power = absolute value of exponent
while power > 0:
    if power is odd: result *= factor
    power = floor(power / 2)
    if power > 0: factor *= factor
return result (or its reciprocal-base computation for negative exponents)
```

For `e = |exponent| > 0`, `T(e) = T(floor(e/2)) + Θ(1)`, so iteration gives
Θ(log e) time in best, average, and worst cases for a fixed nonzero exponent.
The multiplication count is `floor(log₂ e)` squarings plus
`popcount(e)` result multiplications; auxiliary space is Θ(1). For exponent
zero, time and multiplication count are Θ(1) and zero respectively. The API
uses `double`: results are subject to rounding and underflow; non-finite inputs
and overflow to a non-finite intermediate/result are reported as failure.

### Celebrity identification

The adjacency matrix is row-major and any nonzero byte means “knows.” Candidate
elimination compares one matrix entry per person after the initial candidate;
a final scan verifies that everyone knows the candidate and the candidate
knows nobody else. Diagonal entries are ignored.

```text
candidate = 0
for person from 1 to n - 1:
    if candidate knows person: candidate = person
for each other person:
    reject unless person knows candidate and candidate does not know person
return candidate
```

There is no recursive recurrence. For `n > 0`, elimination takes `n - 1`
relationship checks and verification takes at most `2(n - 1)`, yielding Θ(n)
best, average, and worst time after the matrix is available, with Θ(1)
auxiliary space. The reported comparison count is at most `3n - 3`.

### Classical matrix multiplication

Matrices are row-major `double` arrays. For `m × k` and `k × p` inputs, each
output cell is the dot product of one row and one column. Compatible rectangular
dimensions are supported; a zero inner dimension produces a zero-filled
result.

```text
for each row i and output column j:
    result[i][j] = 0
    for k from 0 to inner_dimension - 1:
        result[i][j] += left[i][k] * right[k][j]
```

For positive dimensions, `T(m,k,p) = Θ(mkp)` (best, average, and worst). The
operation count is `mkp` multiplications and `mp(k - 1)` additions when `k`
is positive. It uses Θ(1) working space beyond the Θ(mp) returned matrix.
Dimension arithmetic, allocation sizes, and non-finite arithmetic results are
checked; the output remains empty on failure.

### Strassen matrix multiplication

Strassen requires two equally sized square matrices. The implementation uses
the standard seven half-size products, reuses two input-sum buffers at each
recursive level, and switches to classical multiplication at size 32 or
smaller. Sizes above the cutoff that are not powers of two are padded with
zeroes to the next power of two and cropped back to the requested dimensions.

```text
if size <= 32: use classical multiplication
form the seven Strassen products from quadrant sums/differences
combine the seven products into four result quadrants
```

With the classical base case, `T(n) = 7T(n/2) + Θ(n²)`. Since
`log₂ 7 > 2`, the Master Theorem gives Θ(n^log₂7), approximately Θ(n^2.807),
for best, average, and worst cases. For padded size `p`, the actual bound is
Θ(p^log₂7); since `n ≤ p < 2n`, this is Θ(n^log₂7). Arithmetic additions and
multiplications are also Θ(n^log₂7). The recursive workspace and padding use
Θ(n²) auxiliary space, in addition to the Θ(n²) returned matrix.

Matrix operations use `double` and reject non-finite inputs or any non-finite
intermediate. Floating-point rounding means Strassen and classical results may
differ slightly; callers comparing results should use a tolerance. For small
inputs Strassen deliberately uses the classical method, so it is not faster at
every size.

## Dynamic programming algorithms

The Milestone 4 public API is in
`include/genome/algorithms/dynamic_programming.h`; implementations are in
`src/algorithms/dynamic_programming.c`.

### 0/1 knapsack

Weights and capacity are nonnegative `size_t` values; values and the maximum
result use nonnegative `int64_t` values. Inputs are borrowed and unchanged.
`dp_knapsack_01` is the value-only compatibility API;
`dp_knapsack_01_with_selection` also returns selected indices in ascending
order. The caller owns that index array and must release it with `free`; an
empty selection is represented by a null pointer and a zero count.

Both functions use a one-dimensional tabulation array indexed by capacity.
For each item, capacities are visited in descending order so that the item
can contribute at most once. The selection API additionally stores one
decision byte for each item/capacity state. This uses O(nC) additional bytes
for `n` items and capacity `C`, while retaining O(C) value-table space. During
reconstruction, a strictly better candidate selects the current item; ties
keep the previously computed solution and skip that item, making results
deterministic. Zero-weight items are supported, including when capacity is
zero. Empty input succeeds with a maximum value of zero and no selected
indices. Negative values, invalid pointers, allocation-size overflow,
allocation failure, and an unrepresentable result are reported as failure.
Outputs are reset on failure, and partial allocations are freed.

```text
best[0..capacity] = 0
take[0..item_count-1][0..capacity] = false
for each item i:
    for c from capacity down to weight[i]:
        candidate = best[c - weight[i]] + value[i]
        if candidate > best[c]:
            best[c] = candidate
            take[i][c] = true
reconstruct selected indices by walking items backward from capacity
return best[capacity]
```

For `n` items and capacity `C`, both APIs take O(n(C + 1)) time. The
value-only API uses O(C + 1) auxiliary space; reconstruction uses
O(n(C + 1)) additional decision storage and O(C + 1) value-table space,
instead of an O(n(C + 1)) two-dimensional value table.

### Longest common subsequence

`dp_longest_common_subsequence` accepts two NUL-terminated byte strings;
`dp_longest_common_subsequence_n` accepts the same inputs with explicit byte
lengths. The `_n` inputs must provide the specified number of bytes and cannot
contain embedded NUL bytes. These APIs work for DNA strings as well as other
byte strings. On success, the caller owns the reconstructed, NUL-terminated
LCS string and releases it with `free`; the reported length excludes the
terminator. Free a previous LCS result before reusing its output variable;
these functions reset the output pointer but do not free an earlier result.

Let `L[i][j]` be the LCS length of the first `i` bytes of the first input and
the first `j` bytes of the second input. The base cases are `L[0][j] = 0` and
`L[i][0] = 0`. If the last bytes match, an optimal subsequence includes that
byte and `L[i][j] = L[i - 1][j - 1] + 1`. Otherwise, an optimal subsequence
must exclude one of those last bytes, so `L[i][j]` is the larger of
`L[i - 1][j]` and `L[i][j - 1]`.

The recurrence has optimal substructure because an LCS of matching prefixes
is built from an LCS of smaller prefixes. It has overlapping subproblems
because the same prefix pairs are needed by multiple later states. Tabulation
fills a row-major `(m + 1) × (n + 1)` table from the base cases. Reconstruction
walks backward from `L[m][n]`: matching bytes are copied into the output from
right to left and both indices move diagonally; otherwise the index toward
the larger neighboring table value moves. A tie moves up, making the result
deterministic when several LCS strings are valid.

```text
L[0..m][0..n] = 0
for i from 1 to m:
    for j from 1 to n:
        if first[i - 1] == second[j - 1]:
            L[i][j] = L[i - 1][j - 1] + 1
        else:
            L[i][j] = max(L[i - 1][j], L[i][j - 1])

i = m; j = n; output_position = L[m][n]
while i > 0 and j > 0:
    if first[i - 1] == second[j - 1]:
        output[--output_position] = first[i - 1]
        i--; j--
    else if L[i - 1][j] >= L[i][j - 1]:
        i--
    else:
        j--
return L[m][n] and output
```

For input lengths `m` and `n`, tabulation and reconstruction take O(mn)
time. The DP table uses O(mn) space; the returned string uses an additional
O(min(m, n)) bytes, so total space remains O(mn). The implementation checks
dimension and allocation-size arithmetic before allocating the table.

### Matrix Chain Multiplication

`dp_matrix_chain_multiplication` accepts `matrix_count + 1` positive
dimensions. It returns the minimum scalar multiplication cost as `uint64_t`
and an allocated parenthesization such as `((A1 x A2) x A3)`. The caller owns
the string and must free it; free any previous result before reusing its output
variable. Dimensions and outputs are validated, and outputs remain empty on
failure.

Let `C[i][j]` be the minimum cost of multiplying matrices `Ai` through `Aj`,
and `S[i][j]` be the split index for that optimal order. A single matrix has
cost zero: `C[i][i] = 0`. For `i < j`, choose a split `k` and compute
`C[i][k] + C[k + 1][j] + d[i - 1] * d[k] * d[j]`, using one-based matrix
indices and zero-based dimension-array indices. The minimum candidate is
stored with its split. If costs tie, the lowest split is retained, making the
parenthesization deterministic.

Optimal substructure holds because each side of an optimal split must itself
use a minimum-cost order; otherwise that side could be replaced to reduce the
total cost. Overlapping subproblems arise because candidate splits reuse
costs for the same matrix intervals. The bottom-up table computes intervals
in increasing chain length. Reconstruction recursively follows `S` from the
full chain to emit matrix labels and parentheses.

```text
for i from 1 to n:
    C[i][i] = 0
for length from 2 to n:
    for i from 1 to n - length + 1:
        j = i + length - 1
        C[i][j] = minimum over k in [i, j):
            C[i][k] + C[k + 1][j] + d[i - 1] * d[k] * d[j]
            save the minimizing k in S[i][j]
reconstruct(i, j):
    if i == j: emit "A" followed by i
    otherwise: emit "(" + reconstruct(i, S[i][j])
               + " x " + reconstruct(S[i][j] + 1, j) + ")"
```

For `n` matrices, the interval table takes O(n³) time and O(n²) auxiliary
space. The cost and split tables and all intermediate multiplication and
addition operations are checked for overflow. An overflowing candidate is
ignored rather than allowed to wrap; the call fails if the full chain has no
representable `uint64_t` cost. The returned parenthesization is O(n log n)
bytes because matrix labels include their indices, and total space remains
O(n²).

### Optimal Binary Search Tree

`dp_optimal_bst` builds a search tree for strictly increasing, unique integer
keys and their nonnegative `uint64_t` successful-search frequencies. This is
the successful-search-only formulation: it does not accept unsuccessful
search probabilities. Empty input succeeds with a zero cost and a null tree.
The caller owns a returned tree and releases it with
`dp_optimal_bst_destroy`; destroy a prior tree before reusing its output
variable. The result nodes store each key and frequency and expose left and
right children.

A greedy choice of the locally most frequent key is not sufficient: selecting
a root changes the depth, and therefore the cost, of every key in both
subtrees. Dynamic programming tries every possible root for each key interval
and reuses the already computed optimal subtree costs.

Let `C[i][j]` be the minimum cost of a BST containing keys `i` through `j`,
using zero-based indices, and `R[i][j]` its chosen root. An empty interval has
cost zero (`C[i][i - 1] = 0`). For a nonempty interval, define the frequency
sum `W(i,j) = sum(f[i..j])`. The recurrence is

```text
C[i][j] = min over r in [i, j] of:
          C[i][r - 1] + C[r + 1][j] + W(i, j)
```

The base case for a single key follows directly: its two empty subtrees cost
zero, so `C[i][i] = f[i]`. Prefix sums compute each interval weight in O(1):
`P[0] = 0`, `P[k + 1] = P[k] + f[k]`, and
`W(i,j) = P[j + 1] - P[i]`. Roots are examined in increasing key order and
only a strictly lower cost replaces the current best, so ties select the
smallest key in that interval.

```text
P[0] = 0
for i from 0 to n - 1:
    P[i + 1] = P[i] + frequency[i]
for length from 1 to n:
    for i from 0 to n - length:
        j = i + length - 1
        weight = P[j + 1] - P[i]
        for r from i to j:
            left = 0 if r == i else C[i][r - 1]
            right = 0 if r == j else C[r + 1][j]
            candidate = left + right + weight
            retain the lowest-cost candidate and its root r
build(i, j):
    if i > j: return empty child
    r = R[i][j]
    return node(key[r], build(i, r - 1), build(r + 1, j))
```

For any optimal tree, its root is some key `r`; replacing either subtree
with a cheaper BST for the same key interval would lower the total search
cost. Thus both subtrees must be optimal and the recurrence considers the
optimal root. The same intervals occur under multiple candidate roots, so
tabulation avoids recomputing these overlapping subproblems.

The standard implementation takes Θ(n³) time in the best, average, and worst
cases (Θ(1) for empty input) and uses Θ(n²) auxiliary space for cost, root,
and representability tables; the frequency prefix array uses Θ(n) space.
Reconstruction allocates Θ(n) result nodes and uses O(n) recursive stack
space in the worst case. Prefix-frequency sums and candidate cost additions
are checked in `uint64_t`; the call fails cleanly if the total frequencies or
the minimum search cost are not representable. Any overflowing candidate is
excluded rather than allowed to affect the minimum.

| Input case | Time | Auxiliary space |
|---|---:|---:|
| Empty input | Θ(1) | Θ(1) |
| Nonempty input, best / average / worst | Θ(n³) | Θ(n²) |

### Bellman–Ford shortest paths

`dp_bellman_ford` accepts a directed edge list of
`DpBellmanFordEdge { from, to, weight }` records. Its `int64_t` weights may be
positive, zero, or negative. This edge-list type is separate from the
existing `WeightedEdge` because that graph API stores `double` weights, which
cannot represent every `int64_t` value exactly. Unlike Dijkstra, Bellman–Ford
supports negative weights and detects negative-weight cycles reachable from
the selected source. An unreachable negative cycle does not affect the result.

The conceptual DP state `D[k][v]` is the minimum cost of a path from the
source to vertex `v` using at most `k` edges. The source starts at zero;
other vertices are initially unreachable. For each pass, every edge
`(u, v, w)` may improve the destination:

```text
D[0][source] = 0
D[0][v] = unreachable for v != source
for k from 1 to V - 1:
    D[k][v] = D[k - 1][v] for every vertex v
    for each directed edge (u, v, w):
        if D[k - 1][u] is reachable:
            D[k][v] = min(D[k][v], D[k - 1][u] + w)
```

The implementation relaxes edges in place, which can carry improvements
through multiple edges during a pass; this preserves the same minimum while
often converging earlier. It stops after a pass with no changes. Without a
reachable negative cycle, a shortest path can be chosen simple: any repeated
vertex would form a cycle, and removing a nonnegative cycle cannot increase
the path cost. A simple path contains at most `V - 1` edges, so `V - 1`
complete passes are sufficient. After those passes, an additional scan
reports a negative cycle if a reachable edge can still improve a distance.

```text
reachable[source] = true; distance[source] = 0
for every other vertex v:
    reachable[v] = false
    predecessor[v] = none
repeat at most V - 1 times:
    changed = false
    for each edge (u, v, w):
        if u is reachable and distance[u] + w improves v:
            update distance[v] and predecessor[v]
            mark v reachable; changed = true
    if not changed: stop
for each edge (u, v, w):
    if u is reachable and it can still improve v:
        report a reachable negative cycle
```

The relaxation recurrence is optimal substructure: a shortest path ending
with edge `(u,v)` consists of a shortest path to `u` followed by that edge.
The DP states overlap because many candidate paths reuse the same best
distance to a vertex. `DpBellmanFordResult` owns distance, predecessor, and
reachability arrays; initialize it to `{0}`. A successful call replaces and
destroys the prior result, while a failed call preserves it. Unreachable
distances are unspecified and must be checked through `reachable`. The
predecessor of the source and of unreachable vertices is `SIZE_MAX`.

Path reconstruction returns a caller-owned vertex-index array, released with
`free`. It validates the destination and reachability, limits traversal to at
most `V` vertices to reject cyclic predecessor chains, and refuses
reconstruction when a reachable negative cycle was reported.

| Case | Time | Auxiliary space |
|---|---:|---:|
| Best case with early convergence | O(V + E) | O(V) |
| Average case | Depends on graph and edge order; O(VE) upper bound | O(V) |
| Worst case | O(VE) | O(V) |

All distance additions are checked before evaluation. If a candidate path
distance cannot be represented by `int64_t`, the call fails and preserves
the caller's previous result rather than wrapping. Dijkstra remains the
choice for nonnegative weights; Bellman–Ford is required for negative weights
and reachable negative-cycle detection.

### Floyd–Warshall all-pairs shortest paths

`dp_floyd_warshall` accepts a directed, row-major `V × V` matrix through
`DpWeightedDirectedMatrix`. The parallel `edge_present` array distinguishes an
absent edge from a valid zero-weight edge; any nonzero byte means that the
edge exists. Weights may be positive, zero, or negative, and self-loops are
permitted. The API does not use a numeric sentinel for missing edges.

The DP state `D[i][j]` is the minimum known path cost from vertex `i` to
vertex `j`. Each intermediate-vertex stage either retains the previous path
or joins a path from `i` to `k` with one from `k` to `j`:

```text
D[i][j] = weight(i, j) if edge (i, j) exists
D[i][i] = 0
D[i][j] = unreachable otherwise
for k from 0 to V - 1:
    for i from 0 to V - 1:
        for j from 0 to V - 1:
            if D[i][k] and D[k][j] are reachable:
                D[i][j] = min(D[i][j], D[i][k] + D[k][j])
```

In the recurrence, stage `k` permits only vertices `0..k` as intermediates.
An optimal path either omits `k` or splits at `k`; these are the two
subproblems considered. The same endpoint pairs are reused across stages,
giving overlapping subproblems and optimal substructure. The diagonal is
initialized to zero, except that a negative self-loop is retained and treated
as a negative cycle.

The implementation keeps a next-hop table and reachability table alongside
the distances. A negative diagonal after the DP identifies a negative cycle,
including one disconnected from other vertices. For each such cycle vertex
`k`, every pair `(i,j)` with `i` able to reach `k` and `k` able to reach `j`
is marked in `negative_cycle_affected`. The global flag indicates that at
least one negative cycle exists; the per-pair flags identify exactly which
shortest paths are undefined. Paths that are unreachable or affected by a
negative cycle cannot be reconstructed. Other pairs remain usable even when
some unrelated negative cycle exists.

```text
for every pair (i, j):
    initialize distance and next hop from the matrix and diagonal rule
for k from 0 to V - 1:
    for i from 0 to V - 1:
        for j from 0 to V - 1:
            if both subpaths are reachable:
                checked_candidate = distance[i][k] + distance[k][j]
                if (i, j) is unreachable or checked_candidate is smaller:
                    update distance[i][j]
                    next_hop[i][j] = next_hop[i][k]
for every k with distance[k][k] < 0:
    mark all reachable pairs through k as negative-cycle-affected
```

`DpFloydWarshallResult` owns its distance, next-hop, reachability, and
negative-cycle-affected arrays. Initialize it to `{0}` and destroy it with
`dp_floyd_warshall_result_destroy`. A successful call replaces the previous
result; validation, allocation, or arithmetic failure preserves it.
Reconstructed vertex-index paths are caller-owned and released with `free`.
Path reconstruction is bounded by `V` hops to reject cyclic next-hop chains.

The algorithm takes Θ(V³) time in best, average, and worst cases and Θ(V²)
auxiliary space. All matrix-size calculations and signed distance additions
are checked. If an evaluated reachable candidate cannot fit in `int64_t`,
the computation fails rather than wrapping.

Bellman–Ford computes distances from one source in an edge-list graph and
reports only negative cycles reachable from that source. Floyd–Warshall
computes all-pairs distances from an adjacency matrix and detects negative
cycles anywhere in the graph, while marking the affected source/destination
pairs. On graphs without negative cycles, repeated Bellman–Ford runs from
every vertex provide an independent distance comparison for each Floyd row.

### Travelling Salesman Problem with bitmask DP

`dp_tsp_bitmask` solves the directed, asymmetric TSP over a
`DpWeightedDirectedMatrix`. As for Floyd–Warshall, `edge_present` marks
existing edges, so zero and negative finite costs are valid and cannot be
confused with a missing edge. The default maximum is
`DP_TSP_MAX_VERTICES == 18`; it can be set at compile time from 1 through 20.
This upper bound limits exponential state storage. At the maximum configurable
size, the cost and parent tables together require roughly 188 MiB, excluding
the input and result.

The state `dp[mask][last]` is the minimum cost of a path that starts at the
requested vertex, visits exactly the vertices in `mask`, and ends at `last`.
The base state is `dp[1 << start][start] = 0`. For each absent vertex `next`
connected from `last`, the transition is:

```text
dp[mask | (1 << next)][next] =
    min(dp[mask | (1 << next)][next],
        dp[mask][last] + cost[last][next])
```

For the full mask, add a valid return edge from `last` to `start` and choose
the least representable complete-tour cost. A one-vertex tour is represented
as `[start, start]` and requires an explicit self-loop. A graph with no
Hamiltonian cycle returns `DP_TSP_NO_TOUR`; invalid input, allocation failure,
or an overflowing evaluated candidate returns `DP_TSP_ERROR`.

```text
dp[1 << start][start] = 0
for each state (mask, last) containing start and last:
    for each unvisited next with an edge last -> next:
        candidate = checked_add(dp[mask][last], cost[last][next])
        retain candidate and parent if it improves the destination state
best = minimum over last of:
       dp[full_mask][last] + cost[last][start], when the closing edge exists
reconstruct parents backward and close the tour at start
```

Optimal substructure holds because a minimum path to `(mask,last)` must end
with a minimum path to its preceding state plus the final edge. Distinct
visited-set/final-vertex combinations recur across candidate tours, so the
bitmask table reuses overlapping subproblems. The tour contains the requested
start at both ends, visits each vertex once in between, and is caller-owned;
release it with `free`. Initialize `DpTspResult` to `{0}` and destroy it with
`dp_tsp_result_destroy`. A successful tour or valid no-tour result replaces
the previous result; an error preserves it.

There are O(n·2ⁿ) states and at most O(n) transitions per state, giving
O(n²·2ⁿ) time and O(n·2ⁿ) auxiliary space. This is exponential because every
subset of visited vertices must be distinguished; increasing `n` by one
approximately doubles the number of subsets. Checked bit shifts,
state-table dimensions, allocations, and signed cost additions prevent
undefined behavior or wrapped costs. A single overflow encountered while
evaluating a reachable candidate causes an error, even if another tour might
have a representable cost.

## Backtracking and Branch and Bound

The Milestone 5 API is separate from Dynamic Programming in
`include/genome/algorithms/backtracking_algorithms.h`; implementations are in
`src/algorithms/backtracking_algorithms.c`. All input arrays are borrowed.
Owning results must be initialized to `{0}`. A successful call replaces and
destroys an existing result, while a failed call preserves it. Destroy Sum of
Subsets and Knapsack results with their matching destroy functions.

### Sum of Subsets

`backtracking_sum_of_subsets` accepts strictly positive `int64_t` values and
a nonnegative target. The output contains each index selection exactly once.
Equal values at different input positions are distinct items, so different
index subsets are retained even if their values match. Returned index lists
are sorted in ascending original-input order. For target zero, the empty
subset is the sole answer, including for empty input; for empty input and a
positive target, the result contains no subsets.

The search sorts a private copy by value and then original index. Since all
values are positive, branches stop when the next sorted value exceeds the
remaining target. Every selection is reached by choosing increasing positions,
so no selection is duplicated.

```text
search(first, remaining):
    if remaining == 0: emit the current index set; return
    for i from first to n - 1:
        if value[i] > remaining: stop
        add original_index[i] to current set
        search(i + 1, remaining - value[i])
        remove original_index[i]
```

There are O(2ⁿ) candidate subsets in the worst case, plus O(n log n) sorting.
Recursion and current-path storage use O(n) auxiliary space. Returned memory
is proportional to the number and total sizes of all emitted subsets. Positive
inputs make the remaining-target subtraction safe; invalid nonpositive values,
negative targets, allocation failure, and allocation-size overflow fail
without replacing an existing result.

### 0/1 Knapsack by backtracking

`backtracking_knapsack_01` is independent of the existing
`dp_knapsack_01` implementation. It uses the same borrowed `size_t` weights,
nonnegative `int64_t` values, capacity, and comparable maximum-value/index
result. Each level decides whether to omit or include the next item. An
include branch is created only if its weight fits the remaining capacity.
Zero-weight items are ordinary include decisions and are supported at zero
capacity.

```text
search(i, weight, value):
    if i == n:
        keep this feasible selection if its value improves the best
        return
    search(i + 1, weight, value)
    if item_weight[i] <= capacity - weight:
        search(i + 1,
               weight + item_weight[i],
               value + item_value[i])
```

Every feasible subset corresponds to one root-to-leaf path; infeasible include
branches are omitted. Therefore the maximum among visited leaves is optimal.
The worst-case time is O(2ⁿ); recursion and selection paths use O(n) auxiliary
space, excluding the returned indices. A feasible value sum outside `int64_t`
is reported as failure rather than wrapped.

### 0/1 Knapsack by Branch and Bound

`branch_bound_fifo_knapsack_01`,
`branch_bound_lifo_knapsack_01`, and
`branch_bound_least_cost_knapsack_01` solve the same problem independently
using FIFO queue, LIFO stack, and a max-priority heap keyed by the optimistic
value bound. The heap order is equivalent to choosing the least cost when
cost is defined as the negative upper bound. Nodes carry the number of sorted
items already decided, current weight/value, parent node, and the decision
that created them. Parent indices in a node arena reconstruct the selected
input indices without copying a decision array into every node.

Items are sorted by descending value-to-weight ratio using exact rational
comparison (not floating point), then by original index. Positive-value
zero-weight items sort first. For the remaining items after level `l`, let
`S[l]` be their total value, `Z[l]` their zero-weight value, `r[l]` their
largest positive-weight value/weight ratio (or zero when none remains), and
`R` the remaining capacity. The bound is

```text
current_value + min(S[l], Z[l] + ceil(r[l] * R))
```

Any feasible selection's positive-weight value is at most its maximum
remaining density times its total weight, and that weight cannot exceed `R`.
The zero-weight values are added separately. The total-value suffix is also
an upper bound, so the minimum of these two upper bounds remains optimistic.
The implementation calculates the ceiling using overflow-safe integer
arithmetic and saturates at `INT64_MAX`; it does not use a rounded floating
point ratio. Bounds are precomputed from suffix summaries and evaluated in
O(1). A node is pruned only when its bound is strictly less than the best
value, or equal when the best value is representable below `INT64_MAX`.
At `INT64_MAX`, tied nodes are still explored so a feasible overflowing
selection cannot be hidden by saturation.

```text
sort items by decreasing value / weight, preserving original indices
initialize root with value 0 and its optimistic bound
while the queue, stack, or heap is not empty:
    select a node according to the chosen policy
    update the best feasible selection
    prune only if its optimistic bound cannot improve the best
    generate the omit child
    generate the include child only when it fits
    calculate each child's bound and enqueue it if it may improve the best
reconstruct the best selected indices from parent links
```

The bound is admissible, so no branch containing a better solution is pruned.
Different node-selection policies change exploration order and usually the
amount of pruning, not the optimum. Equal ratios use original input order;
equal heap priorities use bound, current value, then insertion order for
deterministic behavior. Each search can generate O(2ⁿ) nodes; the O(1) bound
makes worst-case time O(n log n + 2ⁿ), including item sorting. The node arena
and open set together require O(2ⁿ) worst-case auxiliary space. In favorable
instances, pruning can substantially reduce generated nodes.

### Choosing an algorithm

| Technique | Main idea | Typical trade-off |
|---|---|---|
| Backtracking | Explore decision choices and reject infeasible partial states | Simple exhaustive search; may visit most of the state-space tree |
| Branch and Bound | Explore partial solutions and prune using an optimistic objective bound | Preserves exactness while often reducing search; queue, stack, and heap orders affect work |
| Dynamic Programming | Store and reuse results for overlapping subproblems | Predictable table-based work and storage; useful when state dimensions remain manageable |
| Greedy | Commit to a locally best choice without reconsidering it | Often fast; globally optimal only when the problem has the required greedy-choice property |

## Greedy algorithms

The Milestone 3 public API is in
`include/genome/algorithms/greedy_algorithms.h`. Implementations are in
`src/algorithms/greedy_algorithms.c`.

Activity selection sorts a copy of the activities by finish time, then retains
an activity when its start time is at or after the last selected finish time.
Invalid intervals are rejected. Empty input succeeds with a zero selection.

Fractional knapsack sorts items by value-to-weight ratio. Negative values,
negative weights, and negative capacity are rejected; zero-value items are
accepted. Positive-value zero-weight items are included first without
consuming capacity.

Huffman coding builds a binary tree with a min-heap and stores no external
character table. The supported alphabet is up to 256 unique byte symbols.
`greedy_huffman_build` returns a tree owned by the caller, which must release
it with `greedy_huffman_destroy`; build failures release all partial trees.
Encode and decode output buffers are also caller-owned and must be released
with `free`. Failed operations leave output pointers empty and free any partial
buffers. A single-symbol tree uses `0` as the code bit for that symbol;
decoding requires the corresponding zero bits and reproduces the symbol once
per bit.

Dijkstra and Prim use row-major weighted adjacency matrices. `DBL_MAX` denotes
an absent edge. Both require finite matrix entries, and Prim additionally
requires symmetric edge weights. Dijkstra requires non-negative weights and
returns distances and predecessor indices; DBL_MAX is the disconnected
distance. Prim supports finite negative and zero edge weights and returns a
forest that begins at the requested vertex and restarts in each disconnected
component. Kruskal uses an edge-list representation for an undirected graph and
supports any finite edge weight, including negative values.

Dijkstra distances and predecessors and Prim/Kruskal edge arrays are
caller-owned result buffers; release them with their corresponding
`*_destroy` function. Huffman encoded/decoded byte arrays are caller-owned and
released with `free`.

Kruskal edges are sorted in increasing weight order; disjoint-set union uses
path compression and union by rank. Self-loops and edges joining
already-connected components are skipped, so duplicate or reversed duplicate
edges cannot create cycles. The result contains a minimum spanning forest, its
total weight, and a `connected` flag; the flag is false for disconnected and
empty graphs and true for a single vertex. For `V` vertices and `E` input
edges, sorting takes O(E log E), and disjoint-set processing takes
O(E α(V)) amortized; total time is O(E log E) and auxiliary storage is
O(V + E), including the sorted edge copy.

## Flow networks and bipartite matching

The Milestone 6 public interfaces are in
`include/genome/graph/ford_fulkerson.h` and
`include/genome/graph/bipartite_matching.h`; implementations are in
`src/graph/`. Both graph inputs borrow their edge arrays for the duration of a
call. Initialize result structures to `{0}` before first use. Success replaces
and destroys the old result; invalid input or allocation failure leaves it
unchanged. Result arrays are owned by the caller and must be released with
their matching `*_result_destroy` function. Do not copy owning results by
value.

### Ford–Fulkerson with Edmonds–Karp search

`FordFulkersonNetwork` represents directed arcs as
`FordFulkersonEdge { from, to, capacity }`. Capacities are signed `int64_t`
values so the API can reject negative input explicitly. Parallel arcs are
accepted as independent input edges. Source and sink must be valid, distinct
vertices.

The implementation repeatedly finds a shortest (fewest-edge) residual path
using BFS. Each input arc has a forward residual arc and a reverse residual
arc. Augmenting a path subtracts its bottleneck from each forward residual
capacity and adds it to the paired reverse residual capacity. These reverse
arcs allow the algorithm to cancel earlier choices and reroute flow. The
returned `edge_flows` array corresponds to input-edge order. After no
source-to-sink residual path remains, a residual BFS identifies the source
side of a minimum cut. `cut_edges` contains original arcs crossing from that
side to the other, and `cut_capacity` is their total capacity.

```text
flow = 0
initialize each forward residual capacity from its input edge
initialize each reverse residual capacity to 0
while BFS finds a residual source-to-sink path:
    bottleneck = minimum residual capacity on the path
    subtract bottleneck from forward arcs and add it to reverse arcs
    flow += bottleneck
source_side = vertices reachable from source in the final residual graph
cut_edges = original edges leaving source_side
```

Edmonds–Karp terminates because BFS ensures each augmentation uses a shortest
residual path; the standard bound is O(VE²) time. This implementation uses
adjacency lists of paired residual arcs, so auxiliary storage is O(V+E), where
E counts input arcs. It checks size calculations and every residual, flow, and
cut-capacity addition before updating. An unrepresentable result is an error
and does not replace the caller's prior result.

By the max-flow/min-cut theorem, when the algorithm terminates the flow value
equals the capacity of the cut induced by residual reachability. Tests check
this equality, per-edge capacity constraints, and conservation at every
nonterminal vertex.

### Maximum bipartite matching

`BipartiteGraph` stores the left and right partition sizes and a borrowed
`BipartiteEdge { left, right }` list. Every edge must connect an in-range left
vertex to an in-range right vertex. Empty partitions and isolated vertices
are valid. Duplicate edges are redundant: they may appear in adjacency
storage but cannot create extra matched pairs or increase the matching size.

The algorithm processes each left vertex and searches for an augmenting path
with DFS. A free right vertex completes the path. If a visited right vertex is
already matched, the search recursively tries to move its matched left vertex
to another right vertex. The assignments are changed only when a complete
augmenting route has been found.

```text
for each left vertex u:
    clear visited-right markers
    augment(u):
        for each adjacent right vertex v:
            if v was visited: continue
            mark v visited
            if v is unmatched or augment(the left vertex matched to v):
                match u with v
                return success
        return failure
```

After processing all left vertices, no augmenting path remains, so the
matching is maximum. The result contains pair records and both direction
maps; `SIZE_MAX` marks unmatched entries. Each vertex occurs in at most one
pair, and all returned pairs correspond to input edges.

The basic algorithm performs O(VE) work in the worst case, with O(V+E)
adjacency and matching storage. The DFS recursion depth is bounded by the
number of left vertices. Allocation-size overflow, invalid endpoints, NULL
input, and allocation failures are rejected while preserving any previous
result. The test suite also compares deterministic small graphs against a
test-only brute-force matching oracle.

## String matching

The public byte-string API is in
`include/genome/string_matching/string_matching.h`, with all four
implementations in `src/string_matching/string_matching.c`. Each call accepts
text and pattern pointers with explicit `size_t` lengths. The inputs are
borrowed and can contain any byte value, including embedded NULs; matching is
not restricted to DNA bases. Pointers must be non-NULL even for zero-length
inputs. An empty pattern matches at all `text_length + 1` valid positions.

`StringMatchResult` contains an allocated positions array, its match count,
and allocated capacity. Positions are zero-based and emitted in increasing
order; overlapping occurrences are included. Initialize the result with
`string_match_result_init` before use and release it with
`string_match_result_destroy`. A successful call replaces the previous result;
invalid inputs, size overflow, or allocation failure preserve it. Partial
arrays and KMP's LPS array are freed on failure.

### Naive matching

The naive algorithm tries each possible starting position and compares bytes
until the pattern differs or is exhausted. It advances by one after a match,
which naturally preserves overlapping matches.

```text
for each valid start position:
    compare text[start..start + pattern_length) with pattern
    if every byte matches: append start
```

For text length `n` and pattern length `m`, worst-case time is O(nm). Its
working space is O(1), excluding the returned positions array.

### Rabin–Karp matching

Rabin–Karp uses a base-256 polynomial hash modulo 1,000,000,007. The initial
pattern and text-window hashes are computed byte by byte. To move the window,
the leading byte's weighted contribution is removed, the remaining hash is
shifted by one byte, and the incoming byte is added. Hash equality is only a
candidate match: `memcmp` verifies the actual window, so collisions never
produce false matches. The bounded byte/base and modulus operations fit in
`uint64_t`.

```text
hash pattern and first window
for each text window:
    if hashes match and bytes match: append its start
    remove the outgoing byte and add the incoming byte to the rolling hash
```

Hash initialization and window movement take O(n + m) expected time. Character
verification after collisions can make worst-case time O(nm). Working space is
O(1), excluding returned matches.

### Knuth–Morris–Pratt

KMP builds an LPS (longest proper prefix that is also a suffix) array for the
pattern. `lps[i]` records the longest proper prefix ending at `i` that is also
a suffix of the pattern prefix ending at `i`. On a mismatch, this value
selects the next candidate pattern position without revisiting text bytes.
After a match, `lps[m - 1]` is used instead of resetting the pattern position
to zero, preserving overlapping matches.

```text
build LPS from pattern prefixes
text_index = 0; pattern_index = 0
while text_index < n:
    if bytes match: advance both indices
    if full pattern matched:
        append text_index - m
        pattern_index = lps[m - 1]
    else if mismatch and pattern_index > 0:
        pattern_index = lps[pattern_index - 1]
    else:
        advance text_index
```

LPS construction takes O(m), and the search takes O(n), for O(n + m) total
time. The LPS array uses O(m) auxiliary space, excluding returned matches.

### Finite-automaton matching

The finite-automaton matcher builds transitions for all 256 byte values and
states zero through the full pattern length. LPS-derived fallback states let
the table be built in O(256m) time; scanning the text takes O(n). Its table
and LPS storage use O(256m) auxiliary space. This includes arbitrary bytes,
and transitions after a full match preserve overlaps.

| Algorithm | Best / expected time | Worst time | Auxiliary space |
|---|---:|---:|---:|
| Naive | O(nm) | O(nm) | O(1), excluding matches |
| Rabin–Karp | O(n + m) expected | O(nm) | O(1), excluding matches |
| KMP | O(n + m) | O(n + m) | O(m), excluding matches |
| Finite automaton | O(256m + n) | O(256m + n) | O(256m), excluding matches |

## Milestone 8 integration

The reusable integration modules are `include/genome/sequence/sequence_statistics.h`,
`include/genome/benchmark/benchmark.h`, `include/genome/export/result_export.h`,
and `include/genome/cli/cli.h`, implemented under the matching `src/`
directories. The C17 command-line entry point in `src/main.c` supports the
interactive menu, `--demo normal|shock|stress`, `--benchmark`, and the legacy
`genome_tool <input.fasta>` form.

### Sequence statistics

`sequence_statistics_compute` accepts a byte sequence and explicit length.
It counts upper- and lowercase A/C/G/T separately and counts other bytes as
ambiguous. GC percentage is defined only for a nonempty sequence without
ambiguous bytes; otherwise its defined flag is false. It borrows the input and
uses O(n) time and O(1) auxiliary space.

### Benchmarking

`benchmark_measure` calls preparation before timing, times only the operation,
then verifies and cleans up outside the measured interval. A record is
published only after every repetition succeeds. The API uses C17
`timespec_get(TIME_UTC)`; wall-clock adjustment and timer resolution are
platform-dependent, so timings are observational and must not be used as
correctness assertions.

`benchmark_run_suite` generates deterministic inputs and verifies results
before recording them. Its 42 records cover sequence sorting/searching,
four string matchers, 0/1 Knapsack, Dijkstra, Prim, Kruskal, and
Ford–Fulkerson. String matching runs on small, medium, large, and adversarial
datasets. For each dataset, the exact same text and pattern are supplied to
Naive, Rabin–Karp, KMP, and the finite automaton; deterministic generation and
pattern placement are performed before timing. Result positions are compared
to the Naive reference and each record carries the match count as its
solution-quality output. Knapsack rows also report the verified optimum.

The memory fields distinguish two scopes: `measured_result_bytes` is the
actual allocated capacity of the returned positions buffer; this is measured
from the result structure, including a known zero-byte result when there are
no matches. `estimated_auxiliary_bytes` models matcher workspace (KMP's LPS
table and the automaton's transition/LPS tables; zero for Naive and
Rabin–Karp). The latter is an estimate based on algorithm allocation sizes,
not allocator overhead. These are not process-wide peak memory or RSS
measurements. Standard C17 provides no portable process-memory API, so this
suite reports no such platform-specific metric. Other algorithms show memory
as unavailable rather than infer a value. Where the source API provides
operation counters, the record can include them. Data generation,
preparation, result validation, and cleanup are not included in each
algorithm's timed interval. All temporary inputs and returned algorithm
results are released by the benchmark module.

### CSV and JSON output

The exporters write sequence statistics, benchmark records, and
algorithm-comparison records using streaming `FILE` output. They validate
record pointers and numeric values, escape CSV quoting and JSON strings, and
report open/write/close failures. They own no caller data and do not retain
file handles. Export filenames are selected by the user or by a named demo;
the menu does not silently choose an output directory.

### Interactive menu and demonstrations

The menu is driven by caller-supplied `FILE` streams through
`genome_cli_run`, allowing the test suite to exercise it without subprocess
or terminal dependencies. It parses complete input lines, validates numeric
choices, gives FASTA records and result buffers explicit cleanup, and returns
success/failure to `main`. Menu categories include DNA statistics, FASTA
display, sorting/searching, representative divide-and-conquer and greedy
examples, LCS, a comparison of the Knapsack solvers, maximum flow,
byte-string matching, benchmarks, and export.

`genome_cli_run_demo` implements three fixed scenarios:

1. **Normal operation:** read `data/samples/demo.fasta`, print sequence
   statistics, sort records, compare all four string matchers, and export
   sequence and comparison results.
2. **Constraint shock:** solve the same 0/1 Knapsack instance under two
   capacities and show how the optimum changes.
3. **Large-data stress:** generate a deterministic one-million-byte text and
   512 synthetic sequences, compare four matchers on three patterns, time and
   verify all twelve searches against Naive, sort the sequence collection, and
   export the benchmark records. A fixed LCG seed makes the data reproducible.

The stress demonstration reports an estimated input/sequence payload
(excluding matcher workspace and allocator overhead) and measured times.
These observations are not memory profiling or asymptotic proofs. The
PowerShell runner `scripts/run_demos.ps1` runs all three scenarios from the
repository root. Scenario-generated CSV/JSON files are written in the current
working directory and can be removed or retained by the user.

### Complexity theory and alignment scope

The project distinguishes P, NP, NP-hard, and NP-complete in
[`complexity-reference.md`](complexity-reference.md). P vs NP remains
unresolved. The coursework documentation contains no selected polynomial
reduction, and this project makes no claim to implement one. Standard
binary-encoded 0/1 Knapsack decision is weakly NP-complete; the O(nW) DP is
pseudopolynomial in capacity. The LCS API compares bytes and should not be
presented as biological sequence alignment.

The [complexity reference](complexity-reference.md) consolidates time and
auxiliary-space bounds for all implemented algorithms. The
[presentation script](presentation-script.md) describes the scenarios and
their intended explanation.
