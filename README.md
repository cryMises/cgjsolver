# cgjsolver

Integer-only Gauss-Jordan solver for up to 10 variables. It prints every row operation as it goes, then checks the answer against the original equations.

- Row operations use integers only (no fractions until the final answer).
- Several modes: replay a seed, search for the best run, fewest steps, smallest numbers, or a size-vs-steps trade-off.
- Every solution is verified (`Check: ok`).

## Download

Get a build from the [latest release](https://github.com/cryMises/cgjsolver/releases/latest):

| File | System |
|---|---|
| `gjsolver-linux` | Linux x86-64 |
| `gjsolver-macos` | Apple Silicon Macs |
| `gjsolver.exe` | Windows 64-bit |

**Linux / macOS:**

```
curl -L -o gjsolver https://github.com/cryMises/cgjsolver/releases/latest/download/gjsolver-linux
chmod +x gjsolver
./gjsolver
```

Use `gjsolver-macos` in the URL on a Mac. If macOS blocks the file, run `xattr -d com.apple.quarantine gjsolver`.

**Windows:** run `gjsolver.exe` from Command Prompt. SmartScreen may warn because the file is unsigned: choose More info, then Run anyway.

**Android (Termux), Raspberry Pi, older Intel Macs:** the downloads will not run there. Build from source (below).

## Build from source

```
gcc -O2 -o gjsolver gjsolver.c
```

One line, straight from GitHub:

```
curl -sL https://raw.githubusercontent.com/cryMises/cgjsolver/main/gjsolver.c | gcc -O2 -x c - -o gjsolver && ./gjsolver
```

- Termux: use `clang` instead of `gcc` (`pkg install clang curl`).
- Windows: use 64-bit MinGW (`gcc -O2 -static -o gjsolver.exe gjsolver.c`). MSVC is not supported, and 32-bit MinGW cannot build it, because the code uses `__int128` and `__builtin_mul_overflow`.

## Input

```
n seed
a1 a2 ... an b
...
```

- Line 1: `n` is the number of variables (1 to 10), then the seed or mode.
- Next `n` lines: one equation each, `n` coefficients followed by the constant. Use `0` for a missing variable.
- Integers only. Clear fractions first (x/2 + y = 3 becomes x + 2y = 6).
- The system must be square with a unique solution, otherwise you get `Singular: no unique solution.`

Run the program with no input to see this guide. To avoid retyping, save the input in a file and run `./gjsolver < input.txt`.

## Seed / modes

| Seed | What it does |
|---|---|
| `> 0` | Replay that exact seed |
| `0` | Random seed |
| `-1` | Search 1,000,000 seeds for the best run within the limits, then replay it |
| `-2` | Fewest steps (exhaustive search over pivot orders) |
| `-3` | Smallest numbers (smallest largest entry, then fewest steps) |
| `-4` | Size vs steps table, then the best plan with entries <= 99999 |

## Example

Input (2x + y = 5, x - y = 1, seed 1):

```
2 1
2 1 5
1 -1 1
```

Output:

```
Seed 1
 Start
      2    1 |    5
      1   -1 |    1
 Pivot x1 (R1)
 1. R2 = (2)R2 - (1)R1
      2    1 |    5
      0   -3 |   -3
 2. (1/3)R2
      2    1 |    5
      0   -1 |   -1
 Pivot x2 (R2)
 3. R1 += (1)R2
      2    0 |    4
      0   -1 |   -1
 4. (1/2)R1
      1    0 |    2
      0   -1 |   -1
 Diagonal
      1    0 |    2
      0   -1 |   -1
 5. -1R2

 Solution
   x1 = 2
   x2 = 1

 Steps: 5 | Max entry: 5 | Check: ok
```

## Limits

Set by constants at the top of `gjsolver.c`:

| Constant | Value | Meaning |
|---|---|---|
| `MAX` | 10 | Largest n (above 10, 64-bit numbers overflow) |
| `LIM_STEPS` | n^2 + n | Most steps accepted by `-1` |
| `LIM_PEAK` | 99999 | Largest entry accepted by `-1` and `-4` |
| `TIME_LIMIT` | 20 s | Search time for `-2`, `-3` and `-4` |

## Notes

- Dense systems with n = 9 or 10 rarely fit within n^2 + n steps and entries <= 99999. Sparse systems often do. The "Within limits" line in `-1` shows how many seeds qualified.
- For `-2`, `-3` and `-4`, "Not proven (limit hit)" means the time or node limit stopped the search early, so a better plan may exist. Results at n >= 9 can vary with machine speed because of the time limit.
- Solutions are always checked, so a plan that is not optimal is still correct.
