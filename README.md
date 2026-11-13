# simplex-degeneracy

GNU Octave based repository for academic studies and empirical benchmarks on degeneracy and anti-cycling pivot rules in the Revised Simplex Method, grounded in Robert J. Vanderbei's *Linear Programming: Foundations and Extensions* (Chapters 1–6).

## Features

- **Two-Phase Revised Simplex Engine** ([`src/simplex_solver.m`](src/simplex_solver.m))
  - Extended sparse matrix formulation $\bar{A} = [A \ \ I_m]$ in canonical form ($\max c^T x \text{ s.t. } Ax \le b, \ x \ge 0$)
  - Single per-iteration sparse/dense `lu(B)` factorization reused across Primal (`B \ b`), BTRAN (`B' \ c_B`), FTRAN (`B \ N(:, j)`), and tie-breaking
  - Auxiliary problem initialization with artificial variable $x_0$ when $\min(b) < 0$
- **Injectable Pivot Rules**
  - [`src/standard_rule.m`](src/standard_rule.m): Dantzig's largest-coefficient rule (baseline, susceptible to cycling)
  - [`src/bland_rule.m`](src/bland_rule.m): Bland's smallest-index rule (Theorem 3.3 finite termination, zero linear-algebra overhead)
  - [`src/lexico_rule.m`](src/lexico_rule.m): Matrix Lexicographic / Perturbation rule (Theorem 3.2, lazy column-by-column tie-breaking via cached LU factors)
- **Benchmarks & Theoretical Validation**
  - [`src/mps_parser.m`](src/mps_parser.m): Sparse triplet-based MPS reader converting NETLIB models into Octave `.mat` matrices
  - [`src/klee_minty.m`](src/klee_minty.m): Vectorized generator for the $n$-dimensional Klee-Minty deformed hypercube ($2^n - 1$ iterations)
  - [`tests/test_theoretical_validation.m`](tests/test_theoretical_validation.m): Automated assertions for 6-dictionary cycling, Klee-Minty worst-case complexity, Strong Duality, and Complementary Slackness
  - [`benchmarks/run_benchmarks.m`](benchmarks/run_benchmarks.m): Comparative benchmark suite exporting `.mat` and `.csv` telemetry

## Requirements

- **GNU Octave** 10.2.0+
- **curl** and **gcc** (only required to download and decompress the NETLIB `degen2`/`degen3` benchmark files)
- **make** (optional, for shortcut targets)

## Setup & Data Download (NETLIB)

Raw NETLIB MPS instances are excluded from version control. To download and preprocess `degen2` and `degen3` into `data/mps/` and `data/processed/`, run [`scripts/download_netlib.sh`](scripts/download_netlib.sh):

```bash
bash scripts/download_netlib.sh
# or using Make:
make data
```

## Quick Start

Run the theoretical validation suite and the fast benchmark suite (`Cycling_6Dict`, `KleeMinty_n5`, `NETLIB_degen2`) via [`main.m`](main.m):

```bash
octave --no-gui main.m
```

Or from inside an interactive Octave session:

```octave
octave:1> run main.m
```

### Additional Modes & [Makefile](Makefile) Targets

```bash
octave --no-gui main.m --test    # Theoretical validation only (< 1s)
octave --no-gui main.m --bench   # Fast benchmarks only
octave --no-gui main.m --full    # Validation + full benchmarks (includes NETLIB_degen3)

# Equivalent Makefile targets:
make test
make bench
make bench-full
```

## Project Structure

```text
.
├── main.m                             # Main entrypoint (validation + benchmarks)
├── Makefile                           # Convenience targets (all, test, bench, data, clean)
├── scripts/
│   └── download_netlib.sh             # Downloads & decompresses NETLIB degen2/degen3
├── src/
│   ├── simplex_solver.m               # Two-phase Revised Simplex solver with LU caching
│   ├── standard_rule.m                # Dantzig largest-coefficient pivot rule
│   ├── bland_rule.m                   # Bland smallest-subscript anti-cycling rule
│   ├── lexico_rule.m                  # Matrix Lexicographic (perturbation) pivot rule
│   ├── mps_parser.m                   # Sparse CSC MPS parser to Vanderbei standard form
│   └── klee_minty.m                   # Klee-Minty deformed hypercube generator
├── tests/
│   └── test_theoretical_validation.m  # Cycling, Klee-Minty, and Duality certification suite
└── benchmarks/
    └── run_benchmarks.m               # Comparative benchmark runner & CSV/MAT exporter
```

## Sample Results (`octave --no-gui main.m`)

| Dataset | Dim ($m \times n$) | Rule | Status | Objective ($\zeta^*$) | Iterations (F1+F2) | Degenerate (%) | Dual Gap | CPU (s) |
| :--- | :---: | :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `Cycling_6Dict` | $3 \times 4$ | Standard | `MAX_ITER` | $0.0000$ | $100$ ($0+100$) | $100.0\%$ | $0.00 \times 10^0$ | $0.017$ |
| `Cycling_6Dict` | $3 \times 4$ | Bland | `OPTIMAL` | $1.0000$ | $7$ ($0+7$) | $85.7\%$ | $0.00 \times 10^0$ | $0.001$ |
| `Cycling_6Dict` | $3 \times 4$ | Lexico | `OPTIMAL` | $1.0000$ | $2$ ($0+2$) | $50.0\%$ | $0.00 \times 10^0$ | $0.001$ |
| `KleeMinty_n5` | $5 \times 5$ | Standard | `OPTIMAL` | $10^8$ | $31$ ($0+31$) | $0.0\%$ | $0.00 \times 10^0$ | $0.005$ |
| `KleeMinty_n5` | $5 \times 5$ | Bland | `OPTIMAL` | $10^8$ | $15$ ($0+15$) | $0.0\%$ | $0.00 \times 10^0$ | $0.003$ |
| `KleeMinty_n5` | $5 \times 5$ | Lexico | `OPTIMAL` | $10^8$ | $31$ ($0+31$) | $0.0\%$ | $0.00 \times 10^0$ | $0.006$ |
| `NETLIB_degen2` | $665 \times 534$ | Standard | `OPTIMAL` | $1435.1780$ | $1904$ ($1200+704$) | $69.7\%$ | $2.27 \times 10^{-13}$ | $1.929$ |
| `NETLIB_degen2` | $665 \times 534$ | Bland | `OPTIMAL` | $1435.1780$ | $9133$ ($712+8421$) | $96.3\%$ | $0.00 \times 10^0$ | $8.283$ |
| `NETLIB_degen2` | $665 \times 534$ | Lexico | `OPTIMAL` | $1435.1780$ | $3186$ ($1073+2113$) | $83.9\%$ | $2.27 \times 10^{-13}$ | $11.009$ |
