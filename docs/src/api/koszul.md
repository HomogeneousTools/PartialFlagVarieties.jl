# Koszul Algebra

Pure algebra for solving long exact sequences arising from Koszul
filtrations. This module knows nothing about flag varieties or bundles —
it operates on `Cohomology` objects and produces new `Cohomology` objects.

Two flavours are provided: **dimension-valued** (`BigInt`) and
**symbolic** (`AffineExpr`), for when the long exact sequence does not
uniquely determine all groups.

## Dimension-valued solvers

```@docs
solve_ses_cohomology
solve_koszul_filtration
```

## Symbolic solvers

```@docs
AffineExpr
is_determined
is_zero_expr
symbolic_variable
solve_ses_cohomology_symbolic
solve_koszul_filtration_symbolic
```
