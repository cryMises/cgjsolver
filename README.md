# cgjsolver

Integer-only Gauss-Jordan solver for up to 10 variables. It prints every row operation as it goes, then checks the answer against the original equations.

- Row operations use integers only (no fractions until the final answer).
- Several modes: replay a seed, search for the best run, fewest steps, smallest numbers, or a size-vs-steps trade-off.
- Every solution is verified (`Check: ok`).

**Contents:** [Linux](#linux) | [macOS](#macos) | [Windows](#windows) | [Android (Termux)](#android-termux) | [Input](#input) | [Modes](#seed--modes) | [Example](#example) | [Limits](#limits) | [Notes](#notes) | [Changes](#changes)

Prebuilt files are on the [latest release](https://github.com/cryMises/cgjsolver/releases/latest):

| File | Runs on |
|---|---|
| `gjsolver-linux` | Linux, x86-64 (64-bit Intel/AMD) |
| `gjsolver-macos` | macOS, Apple Silicon (M1 and newer) |
| `gjsolver.exe` | Windows, 64-bit (x64; ARM Windows via emulation) |

A prebuilt file only runs on the system and CPU type it was built for. Anywhere else, build from source.

---

## Linux

### Option 1: prebuilt binary (x86-64 only)

```
curl -L -o gjsolver https://github.com/cryMises/cgjsolver/releases/latest/download/gjsolver-linux
chmod +x gjsolver
./gjsolver
```

On a Raspberry Pi or other ARM Linux machine this fails with `Exec format error`. Use Option 2.

### Option 2: compile from source

1. Install a compiler if you don't have one:

   | Distro | Command |
   |---|---|
   | Debian / Ubuntu | `sudo apt install gcc curl` |
   | Fedora | `sudo dnf install gcc curl` |
   | Arch | `sudo pacman -S gcc curl` |

2. Compile and run:

```
curl -sL https://raw.githubusercontent.com/cryMises/cgjsolver/main/gjsolver.c | gcc -O2 -x c - -o gjsolver && ./gjsolver
```

Or, with the file already downloaded: `gcc -O2 -o gjsolver gjsolver.c`

---

## macOS

### Option 1: prebuilt binary (Apple Silicon only)

```
curl -L -o gjsolver https://github.com/cryMises/cgjsolver/releases/latest/download/gjsolver-macos
chmod +x gjsolver
./gjsolver
```

If macOS refuses to open the file, remove the quarantine flag:

```
xattr -d com.apple.quarantine gjsolver
```

Or open **System Settings > Privacy & Security** and choose **Open Anyway** after the first attempt.

This file will not run on Intel Macs. Use Option 2.

### Option 2: compile from source (Apple Silicon and Intel)

1. Install Apple's command line tools (includes the compiler):

```
xcode-select --install
```

2. Compile and run:

```
curl -sL https://raw.githubusercontent.com/cryMises/cgjsolver/main/gjsolver.c | clang -O2 -x c - -o gjsolver && ./gjsolver
```

Or, with the file already downloaded: `clang -O2 -o gjsolver gjsolver.c`

On macOS `gcc` is an alias for clang, so either name works. Do not use `-static` on macOS; it is not supported.

---

## Windows

### Option 1: prebuilt binary (64-bit)

1. Download `gjsolver.exe` from the [latest release](https://github.com/cryMises/cgjsolver/releases/latest).
2. Open **Command Prompt** or **PowerShell** in the folder where you saved it.
3. Run it:
   - Command Prompt: `gjsolver.exe`
   - PowerShell: `.\gjsolver.exe`

Don't double-click it: the window closes as soon as the program finishes.

Windows SmartScreen may warn because the file is unsigned. Choose **More info**, then **Run anyway**.

### Option 2: compile from source

1. Install a 64-bit MinGW-w64 `gcc`. The simplest way is MSYS2:
   - Install MSYS2 from https://www.msys2.org
   - Open the **MSYS2 UCRT64** terminal and run: `pacman -S mingw-w64-ucrt-x86_64-gcc`
   - Either work inside that terminal, or add `C:\msys64\ucrt64\bin` to your Windows `PATH` to use `gcc` from Command Prompt and PowerShell.

2. Compile and run:

   **Command Prompt:**
   ```
   curl -sL https://raw.githubusercontent.com/cryMises/cgjsolver/main/gjsolver.c | gcc -O2 -static -x c - -o gjsolver.exe && gjsolver.exe
   ```

   **PowerShell** (`curl` is an alias there, so use `curl.exe`, and avoid the pipe):
   ```
   curl.exe -sLO https://raw.githubusercontent.com/cryMises/cgjsolver/main/gjsolver.c; gcc -O2 -static -o gjsolver.exe gjsolver.c; .\gjsolver.exe
   ```

   **File already downloaded:** `gcc -O2 -static -o gjsolver.exe gjsolver.c`

`-static` builds the runtime into the exe so it does not need extra DLLs on other machines.

**Not supported on Windows:** MSVC (Visual Studio's compiler) and 32-bit MinGW. The code uses `__int128` and `__builtin_mul_overflow`, which they lack.

---

## Android (Termux)

The release downloads do not run on Android (the CPU is ARM). Compile instead:

```
pkg install clang curl
curl -sL https://raw.githubusercontent.com/cryMises/cgjsolver/main/gjsolver.c | clang -O2 -x c - -o gjsolver && ./gjsolver
```

---

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

Run the program with no input to see this guide. To avoid retyping, save the input in a file:

| Shell | Command |
|---|---|
| Linux, macOS, Termux, Command Prompt | `./gjsolver < input.txt` (Windows: `gjsolver.exe < input.txt`) |
| PowerShell | `Get-Content input.txt \| .\gjsolver.exe` |

## Seed / modes

| Seed | What it does |
|---|---|
| `> 0` | Replay that exact seed |
| `0` | Random seed |
| `-1` | Search 1,000,000 seeds for the best run within the limits, then replay it. The limits are the bottom 2.5% of steps and of largest entry over the seeds searched |
| `-2` | Fewest steps (exhaustive search over pivot orders) |
| `-3` | Smallest numbers (smallest largest entry, then fewest steps) |
| `-4` | Size vs steps table, then the fewest-step plan with entries <= max(100, 10^(2n - 6.5)) |

### How `-1` picks a seed

It prints a steps table and a largest-entry table for the seeds it searched. Each is split into five groups (bottom 2.5%, bottom 16%, middle 50%, top 84%, top 97.5%) using the mean and the mean plus or minus 1 and 2 standard deviations of that run's own seeds. The entry table works on log10 of the largest entry.

The limits are the **bottom 2.5%** rows of those two tables:

1. If some seeds are within both limits, it picks the one with the fewest steps, then the smallest entry.
2. If none are, it picks the seed **closest to both limits**. Each overshoot is measured in standard deviations, and the combined distance is minimised. The `SEEDS` section shows it as `Closest to both`.

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
| `LIM_STEPS`, `LIM_PEAK` in `-1` | bottom 2.5% rows of the run's own tables | Set after the seeds are run, not by a constant |
| `LIM_PEAK` in `-4` | max(100, 10^(2n - 6.5)) | Largest entry accepted by `-4` (from `PEAK_FLOOR`, `PEAK_SLOPE`, `PEAK_SHIFT`) |
| `BIG` | 2,000,000,000 | `-4` drops any plan whose entries pass this |
| `TIME_LIMIT` | 20 s | Search time for `-2`, `-3` and `-4` |

## Notes

- Because the `-1` limits are the bottom 2.5% rows, few or no seeds may be within both. The `LIMITS` section in `-1` shows how many qualified, and when none did, the closest seed is used.
- At n = 10 most seeds overflow 64 bits and are skipped (the search line says how many), so the `-1` tables come from the rest.
- `-4` has no seed statistics, so it keeps the fixed entry limit above. Its "proven best" holds only within the 2e9 cap, so `-2`, which has no cap, can occasionally find fewer steps.
- For `-2`, `-3` and `-4`, "Not proven (limit hit)" means the time or node limit stopped the search early, so a better plan may exist. Results at n >= 9 can vary with machine speed because of the time limit.
- Solutions are always checked, so a plan that is not optimal is still correct.
- Only run one-line install commands on code you trust: they compile and run whatever the link serves.

## Changes

- **`-1` limits:** the steps and entry limits are now the bottom 2.5% rows of the run's own tables, replacing the fixed formulas. The usage text printed by the program changed to match.
- **`-1` fallback:** if no seed is within both limits, it picks the seed closest to both, instead of just the one with the fewest steps.
- **Entry table:** seed counts per row are now taken against the integer threshold printed in that row, so the table and the limits always agree. Thresholds are rounded to the nearest integer.
- **Overflow:** the top 97.5% entry row at n = 10 used to print a garbage value (`-9223372036854775808`) when the threshold passed 64 bits. It is now clamped to `9.22e+18`.
- **Unchanged:** `-2`, `-3`, `-4` and fixed seeds give the same output as before.
- **Docs:** the `Limits` and `Notes` sections described older fixed limits (n^2 + n steps, entries <= 99999), and now match the program.
