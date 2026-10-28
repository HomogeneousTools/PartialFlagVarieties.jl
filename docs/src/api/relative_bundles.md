# Relative Grassmann bundles

`GrassmannBundle(k, E)` parametrizes subspaces of dimension `k` in the fibres
of `E`. The universal sequence is `0 → S → π*E → Q → 0`.
The projective bundle of lines in `E` is `GrassmannBundle(1, E)`, whose
`O(1)` pushes forward to `dual(E)`.
`GrassmannBundle(rank(E), E)` is the identity projection onto the base.

`GrassmannBundle` represents the space; `RelativeGrassmannianBundle`
represents a vector bundle on that space.

Use `pullback`, `O`, and the universal bundle constructors to build these
vector bundles.

```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(2);

julia> E = direct_sum(line_bundle(X, -1), 3 * structure_sheaf(X));

julia> G = GrassmannBundle(2, E);

julia> pushforward(G, dual(universal_subbundle(G)))[0] == dual(E)
true

julia> pushforward(G, O(G, 1))[0] == exterior_power(dual(E), 2)
true
```

Relative Borel--Weil--Bott produces Schur functors of `dual(E)` on the base.
The full general linear group character is retained; determinant twists are
essential for negative weights and for compatibility with twisting `E` by a
line bundle. `schur_functor` supports direct sums and integral dominant weights.

For a direct sum of relative Levi bundles each summand has only one nonzero
higher direct image. Leray therefore degenerates summand by summand.
Bundles carrying a filtration, including the total cotangent bundle, can have
nonzero differentials. Their cohomology remains symbolic when necessary;
Euler characteristics are exact. The universal extension is never identified
with the direct sum of its graded pieces.

The bundle `E` defining the Grassmann bundle may also be filtered.
In that case its Schur functors carry the induced filtration on the base.
For example, `GrassmannBundle(2, universal_quotient_bundle(quadric(3)))`
retains the extension in the rank-four quotient bundle on the quadric.

```@docs
GrassmannBundle
RelativeGrassmannianBundle
pushforward(::GrassmannBundle, ::RelativeGrassmannianBundle)
```
