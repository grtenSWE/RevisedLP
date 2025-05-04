# RevisedLP

**RevisedLP** is a MATLAB implementation of the revised simplex algorithm for solving linear programming problems of the form:

> minimize \( c^T x \) subject to \( Ax = b, x => 0 \)

This project includes a full two-phase simplex method, robust handling of infeasible and unbounded problems, and Bland’s Rule to prevent cycling.

---

## Features

- Two-phase simplex method with artificial variables (Phase I)
- Revised simplex implementation for efficient matrix operations
- Bland’s anti-cycling rule for numerical stability
- Clean separation of solver logic and test harnesses
- Automatic detection of optimality, infeasibility, and unboundedness
- Extensive test cases including:
  - Optimal BFS
  - Infeasible problems
  - Unbounded problems
  - Degenerate bases
  - Large-scale LPs

---

## File Structure

| File | Description |
|------|-------------|
| `fullsimplex.m` | Main driver: two-phase revised simplex algorithm |
| `revisedsimplex.m` | Phase I & II solver using revised simplex |
| `revisedfindenter.m` | Selects entering variable with Bland’s Rule |
| `revisedfindleave.m` | Performs minimum-ratio test for leaving variable |
| `bupdate.m` | Updates basis matrix and index vectors after pivot |
| `simplex.m` | Original simplex method (1-phase) |
| `findenter.m`, `findleave.m` | Helpers for the original simplex method |
| `simplex_tests.m` | Test suite for original simplex solver |
| `benchmark_simplex.m` | Script to benchmark naive vs optimized solvers |

---

## How to Run

1. Open MATLAB and navigate to the project directory.
2. Run `simplex_tests` to verify correctness of the original simplex solver.
3. Run `benchmark_simplex` to compare performance between naive and optimized versions.
4. Run any test cases in `fullsimplextests` to validate the two-phase solver.

---

## Example Usage

```matlab
A = [1 1; 2 1];
b = [4; 5];
c = [1; 2];
m = 2; n = 2;
[z, x, pi, indices, flag] = fullsimplex(A, b, c, m, n);
```

---

## Performance

This implementation achieves a 2×–10× speed-up over a naive simplex solver using:
- MATLAB’s backslash operator for basis solves (`B\v` instead of `inv(B)*v`)
- Vectorized reduced-cost and ratio calculations
- Logical masking and pre-allocation

---
