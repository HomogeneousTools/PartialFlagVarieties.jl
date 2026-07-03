# Koszul Algebra

Pure algebra for solving long exact sequences arising from Koszul
filtrations. This module knows nothing about flag varieties or bundles —
it operates on `Cohomology` objects and produces new `Cohomology` objects.

Two flavours are provided: **dimension-valued** (`BigInt`) and
**symbolic** (`AffineExpr`), for when the long exact sequence does not
uniquely determine all groups.

## Algorithm overview

A **short exact sequence** of sheaves ``0 \to \mathcal{A} \to \mathcal{B} \to \mathcal{C} \to 0`` induces a
long exact sequence in cohomology:

```math
\cdots \to \mathrm{H}^i(X, \mathcal{A}) \to \mathrm{H}^i(X, \mathcal{B}) \to \mathrm{H}^i(X, \mathcal{C}) \xrightarrow{\delta_i} \mathrm{H}^{i+1}(X, \mathcal{A}) \to \cdots
```

Given the dimensions ``a_i = \mathrm{h}^i(X, \mathcal{A})``, ``b_i = \mathrm{h}^i(X, \mathcal{B})``, ``c_i = \mathrm{h}^i(X, \mathcal{C})``
(any two of the three known), the connecting map ``\delta_i`` has rank
``r_i`` satisfying the rank-nullity constraint ``0 \le r_i \le \min(c_i, a_{i+1})``.
When ``b_i`` is the unknown, the relation

```math
b_i = c_i - r_i + a_i - r_{i-1}
```

determines ``b_i`` exactly if and only if each ``r_i`` is forced by the bounds.

**Dimension-valued solver.** When the bounds force all ranks uniquely (which
is the typical case for cohomology on ``\mathrm{G}/\mathrm{P}``), `solve_ses_cohomology` returns
exact `BigInt` dimensions.

**Symbolic solver.** When the bounds leave genuine ambiguity, each unknown rank
is represented as a fresh [`AffineExpr`](@ref) variable ``x_i``. All entries of
the solved cohomology are affine expressions ``c + \sum_i k_i x_i``. Subsequent
constraints — Serre duality, Euler characteristic, or Akizuki–Nakano vanishing
— can equate these expressions and resolve the variables; see
[`is_determined`](@ref).

**Koszul filtration.** The Koszul complex of a rank-``r`` bundle is chopped
into ``r`` nested short exact sequences. `solve_koszul_filtration` (resp.
`solve_koszul_filtration_symbolic`) iterates the SES solver through the
filtration, accumulating the final cohomology.

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

## Worked examples

### Determined short exact sequence

When one of the three sheaves in a short exact sequence
``0 \to \mathcal{A} \to \mathcal{B} \to \mathcal{C} \to 0`` has vanishing higher cohomology,
the connecting maps ``\delta_i`` are forced to be zero, and the
remaining cohomology groups are uniquely determined.

```jldoctest
julia> using PartialFlagVarieties

julia> A = Cohomology{BigInt}(BigInt[3, 0, 0], 2);  # H⁰=3, H¹=H²=0

julia> B = Cohomology{BigInt}(BigInt[7, 2, 0], 2);  # H⁰=7, H¹=2, H²=0

julia> C, exact = solve_ses_cohomology(A, B)
(H⁰ = 4
H¹ = 2, true)
```

Because ``a_1 = a_2 = 0``, every connecting map ``\delta_i \colon C^i \to A^{i+1}``
has rank 0, so ``c_i = b_i - a_i`` and the result is exact: the second
return value is `true`.

### Ambiguous short exact sequence

When the bounds ``0 \le r_i \le \min(c_i, a_{i+1})`` do not force every
connecting-map rank, the symbolic solver introduces fresh variables.

```jldoctest
julia> using PartialFlagVarieties

julia> A = Cohomology{BigInt}(BigInt[3, 2, 0], 2);  # H⁰=3, H¹=2, H²=0

julia> B = Cohomology{BigInt}(BigInt[5, 5, 1], 2);  # H⁰=5, H¹=5, H²=1

julia> C, exact = solve_ses_cohomology(A, B)
(H⁰ = 2
H¹ = 3
H² = 1, false)

julia> var_counter = Ref(0);

julia> C_sym = solve_ses_cohomology_symbolic(A, B, var_counter)
H⁰ = 2 + x_0
H¹ = 3 + x_0
H² = 1

julia> is_determined(C_sym[0])
false

julia> is_determined(C_sym[2])
true
```

The dimension-valued solver picks *some* valid solution (here
``r_0 = 0``, giving ``c_0 = 2``), but returns `exact = false` to
indicate the answer is not unique. The symbolic solver instead
introduces a variable ``x_0 = \mathrm{rk}(\delta_0)`` and expresses
each ``c_i`` as an affine function of ``x_0``. The variable satisfies
``0 \le x_0 \le 2``, so the true
``\mathrm{h}^0(X, \mathcal{C})`` lies in ``\{2, 3, 4\}``.

### Geometric example: hypersurface (fully determined)

For a *hypersurface* zero locus the Koszul complex has a single short
exact sequence ``0 \to \mathcal{L}^\vee \to \mathcal{O}_X \to \mathcal{O}_Z \to 0``
in which both ambient terms are known from Borel–Weil–Bott, so the solver
determines the result with no rank ambiguity. The quintic threefold
``Z \subset \mathbb{P}^4`` is the classic example:

```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> h = hodge_numbers(Z);

julia> h[2, 3]   # h^{1,2}(Z)
101

julia> all(PartialFlagVarieties.is_determined, h)
true
```

All Hodge numbers are `BigInt` (wrapped in trivial `AffineExpr`): the
Koszul solver encountered no ambiguity.

### Geometric example: higher-rank bundle (symbolic Hodge numbers)

For a zero locus of a *higher-rank* bundle the Koszul filtration has
multiple short exact sequences whose connecting-map ranks need not be
forced. Consider the fourfold ``Z`` cut out by two copies of
``\mathrm{Sym}^2 \mathcal{S}^\vee`` on ``\mathrm{Gr}(2,7)``:

```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 7);

julia> S = universal_subbundle(X);

julia> E = direct_sum(symmetric_power(dual(S), 2), symmetric_power(dual(S), 2));

julia> Z = zero_locus(E);  # dim 4

julia> h = hodge_numbers(Z);

julia> h[2, 2]   # h^{1,1}(Z)
8 + x_0

julia> h[2, 3]   # h^{1,2}(Z)
x_0

julia> h[3, 3]   # h^{2,2}(Z)
30 + 2 * x_0

julia> all(PartialFlagVarieties.is_determined, h)
false
```

Here ``x_0`` is the unknown rank of a connecting homomorphism in one of
the Koszul short exact sequences. Serre duality forces ``\mathrm{h}^{p,q} = \mathrm{h}^{n-p,n-q}``,
and the Euler characteristic ``\chi(\mathcal{O}_Z) = 1`` imposes a further linear
relation, but one degree of freedom remains: Grothendieck–Lefschetz pins
down ``\mathrm{h}^{1,1}`` for complete intersections of ample divisors, but ``Z`` is
the zero locus of a higher-rank bundle that does not split into ample line
bundles, so no Lefschetz-type argument applies.
The true value (computable by other means) is ``\mathrm{h}^{1,1}(Z) = 8``,
corresponding to ``x_0 = 0``.

## Internals

Entry-based long exact sequence solvers used by the restriction machinery.

```@docs
PartialFlagVarieties.les_cokernel
PartialFlagVarieties.les_kernel
```
