# Common Workflows

This page organizes the package by task rather than by module name.

## Construct an arbitrary partial flag variety

Use a named constructor when one exists:

```julia
X = Gr(2, 5)
Y = quadric(5)
Z = cayley_plane()
```

Otherwise use [`partial_flag_variety`](@ref):

```julia
X = partial_flag_variety("D5", [2, 5])
marked_dynkin_diagram(X)
```

The quickest structural invariants are:

```julia
dimension(X)
euler_characteristic(X)
betti_numbers(X)
picard_rank(X)
```

## Compute cohomology of an equivariant bundle

The usual pattern is:

```julia
X = Gr(2, 5)
E = exterior_power(tangent_bundle(X), 2)

H = cohomology(E)
H[0]
chi(H)
```

`cohomology(E)` returns dimensions by default; pass `characters=true` when you
want the full ambient representation in each cohomology group.

## Build standard bundles

The public bundle interface is meant to cover most use cases:

```julia
X = Gr(2, 5)

E = twist(exterior_power(Q(X), 2), 1)
F = direct_sum(E, canonical_bundle(X))
```

Here `Q(X)` is the single-letter shorthand for [`universal_quotient_bundle`](@ref);
the full set (`S`, `Q`, `O`, `T`, `E`) is listed under
[Bundle shorthands](conventions.md#Bundle-shorthands).

For multi-step type-A flags, [`tautological_bundles`](@ref) and
[`universal_subbundles`](@ref) are best thought of as **convenient equivariant
building blocks**, not literal geometric tautological flags.

## Study a zero locus

```julia
X = Gr(2, 5)
E = O(X, 1)
Z = zero_locus(E)

codimension(Z)
normal_bundle(Z)
koszul_terms(Z)
cohomology_on_restriction(Z, T(X))
```

This is the right workflow when you want cohomology of restrictions, Euler
characteristics, Hodge numbers, or Fano / Calabi–Yau tests for ``Z``.

## Interpret symbolic Hodge output

For higher-codimension zero loci, the long exact sequences may leave genuinely
undetermined connecting-map ranks. In that case the package returns
[`AffineExpr`](@ref) entries:

```julia
X = Gr(2, 7)
S = universal_subbundle(X)
E = direct_sum(symmetric_power(dual(S), 2), symmetric_power(dual(S), 2))
Z = zero_locus(E)
H = hodge_numbers_symbolic(Z)
```

Use [`is_determined`](@ref) to test whether an entry is numerical or still
contains a symbolic variable.

## Round-trip external labels

If you work with datasets or published classifications that use ZeroLocus62
labels, the package can encode and decode them:

```julia
label = zerolocus62_label(Gr(2, 4))
X = PartialFlagVariety(label)
```

See [Labels](api/labels.md) for the precise conventions and caveats.
