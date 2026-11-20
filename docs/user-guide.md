# User Guide

## Build and launch (GCC/MinGW)

From the project root in PowerShell:

```powershell
cmake -S . -B build -G "MinGW Makefiles"
cmake --build build
ctest --test-dir build --output-on-failure
.\build\genome_tool.exe
```

If `build` was previously configured with another CMake generator, choose a
new build directory or remove the old build directory yourself before
reconfiguring it with MinGW.

## Navigate the CLI

The main menu numbers its 11 project modules as `[01]` through `[11]`;
`[00]` requests application exit and asks for confirmation. Each module has
its own formatted submenu. Choose `[01]` there to run its existing workflow
or `[00]` to return to the main menu. Invalid selections show an error and
leave the menu available for another choice.

On terminals that support ANSI escape sequences, headings use color and the
screen clears when moving between the main menu and a module. Color is
disabled when `NO_COLOR` is set or the terminal reports `TERM=dumb`.

## Run a DNA alignment

In the main menu, choose **[11] DNA Sequence Alignment**, then select:

- **[01] Needleman-Wunsch** to align both complete sequences.
- **[02] Smith-Waterman** to find their best-scoring matching region.
- **[03] View Scoring Parameters** to see the defaults for the next run.
- **[04] Alignment History** to review alignments from this session.
- **[00] Return to Main Menu** to leave the alignment submenu.

Enter the two original DNA sequences, without adding gaps. Each sequence must
be non-empty and contain only A, C, G, and T; lowercase input is accepted.
Invalid sequences are rejected. The program inserts gaps during traceback.

The CLI asks for a match reward, mismatch penalty, and gap penalty. Press
Enter at a prompt to use the current displayed values (initially `+2`, `-1`,
and `-2`); signed integer values may be entered to change the values for this
and subsequent alignments in the session. Values that would overflow
alignment score arithmetic are reported as an alignment error.

The output includes the score, both aligned strings, a `|` under each match,
counts of matches/mismatches/gaps, and percentage identity (matches divided
by alignment columns). Local alignments also report the matching input ranges
using 1-based inclusive positions. If no positive-scoring local region exists,
the output says so and returns an empty local alignment.

## Examples

### Global alignment

Select `[01] Needleman-Wunsch`, enter `ACG` and `AG`, then press Enter for all
three scoring prompts:

```text
Algorithm       : Needleman-Wunsch
Alignment Type  : Global Alignment
Sequence 1      : ACG
Sequence 2      : AG

Aligned Sequence 1 : ACG
                     | |
Aligned Sequence 2 : A-G

Alignment Score : 2
Matches         : 2
Mismatches      : 0
Gaps            : 1
Identity        : 66.67%
```

The score is `+2 - 2 + 2 = 2`. The complete inputs are represented, and the
gap is generated automatically.

### Local alignment

Select `[02] Smith-Waterman`, enter `TTACGTAA` and `GGACGCCA`, then accept
the defaults:

```text
Algorithm       : Smith-Waterman
Alignment Type  : Local Alignment
Sequence 1      : TTACGTAA
Sequence 2      : GGACGCCA

Aligned Sequence 1 : ACG
                     |||
Aligned Sequence 2 : ACG

Alignment Score : 6
Matches         : 3
Mismatches      : 0
Gaps            : 0
Identity        : 100.00%
Matching region in sequence 1: 3-5 (1-based, inclusive)
Matching region in sequence 2: 3-5 (1-based, inclusive)
```

This local result contains the best-scoring region; the unaligned flanks are
not shown in the aligned strings.

## How the algorithms differ

- **LCS (Longest Common Subsequence)** finds a longest ordered sequence of
  shared characters. It does not penalize mismatches or gaps and does not
  produce a full alignment.
- **Needleman–Wunsch** uses dynamic programming to score an end-to-end
  alignment of both inputs, including substitutions and inserted gaps.
- **Smith–Waterman** uses local resets to zero and traceback from the
  best-scoring cell, returning only the best local region.

Both alignment algorithms use O(mn) time and O(mn) dynamic-programming space
for sequence lengths `m` and `n`.
