# Short Project Presentation Script

## Suggested 4-minute walkthrough

**Opening (20 seconds)**  
“This is the Genome Sequence Search and Alignment Toolkit. It is implemented
in standard C17 and groups the algorithms behind reusable headers and
functions. The command-line menu connects the existing algorithm milestones
to FASTA input, deterministic demonstrations, and benchmarks.”

**Normal operation (45 seconds)**  
Run `scripts\run_demos.ps1` or `genome_tool.exe --demo normal`. Point out the
loaded sample FASTA sequences, base counts and GC percentage, ordering by
length, and the agreement between Naive, Rabin–Karp, KMP, and finite-automaton
matching.

**Constraint shock (35 seconds)**  
Run `genome_tool.exe --demo shock`. The same 0/1 Knapsack instance is solved
with capacity 30 and then capacity 12. Explain that feasibility changes as
the capacity constraint tightens, while the DP implementation still returns
the optimal value.

**Large-data stress (55 seconds)**  
Run `genome_tool.exe --demo stress`. It creates a deterministic 1,000,000-byte
DNA text and 512 synthetic sequences, searches three patterns with all three
string matchers, checks each output against Naive matching, times the
operations, and sorts the sequence collection. The twelve verified search
timings are printed. Explain that elapsed times are observations, not proofs
of asymptotic complexity.

**Complexity theory (60 seconds)**  
“P is the class of decision problems solvable in polynomial time. NP is the
class whose yes answers have certificates verifiable in polynomial time.
NP-hard means every problem in NP has a polynomial-time reduction to the
problem; NP-complete means both NP-hard and in NP. Whether P equals NP is
unknown. For 0/1 Knapsack, the standard decision problem is weakly
NP-complete. The toolkit's O(nW) dynamic program is pseudopolynomial because
W is a numeric capacity, not its encoded bit length. Our coursework
documents do not identify a particular reduction, so I am not claiming to
have implemented or demonstrated one.”

**Close (25 seconds)**  
“The full automated suite protects the earlier milestones. The complexity
reference lists time and auxiliary-space bounds, and the architecture
document explains data ownership and algorithm choices. These results are
algorithm demonstrations; the bytewise LCS function is not a biological
sequence-alignment model.”

## Commands

From the project root after building:

```powershell
cmake -S . -B build
cmake --build build --config Release
ctest --test-dir build -C Release --output-on-failure
.\build\Release\genome_tool.exe
.\build\Release\genome_tool.exe --demo normal
.\build\Release\genome_tool.exe --demo shock
.\build\Release\genome_tool.exe --demo stress
.\build\Release\genome_tool.exe --benchmark
```

For MinGW, the executable is normally `build\genome_tool.exe`. The PowerShell
demo runner selects either executable location automatically.
