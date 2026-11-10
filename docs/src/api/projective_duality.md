# Projective duality

Compute the degree, codimension, and defect of the projective dual of a smooth
partial flag variety or smooth regular zero locus over the complex numbers.
Each invariant has its own function, with an optional polarization:

```jldoctest
julia> using PartialFlagVarieties

julia> X = Gr(2, 5);

julia> projective_dual_degree(X)
5

julia> projective_dual_codimension(X)
3

julia> projective_dual_defect(X)
2
```

The defect is the codimension minus one.
For a hypersurface dual it is zero.
All three functions share a private cache and reuse the Euler characteristics
of the linear sections already computed.

## Choosing the embedding

Omitting `L` uses [`polarization`](@ref) of the source.
An explicit `L` overrides this choice for one call:

```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(2); set_polarization!(X, 3);

julia> projective_dual_degree(X)  # discriminant of ternary cubics
12

julia> projective_dual_degree(X, line_bundle(X, 2))
3

julia> A = projective_space(4); Z = zero_locus(line_bundle(A, 5));

julia> projective_dual_degree(Z)
320

julia> projective_dual_degree(Z, polarization(Z))
320
```

For a zero locus, `L` can be an ambient line bundle or its recorded restriction,
such as `line_bundle(Z, 1)` or `restrict(Z, line_bundle(A, 1))`.
The implementation requires ampleness on the ambient partial flag variety,
which ensures that the restriction is very ample.
It assumes that the defining section is regular and its zero locus is smooth.
It does not test these geometric hypotheses.

The dual is understood in the source's linear span.
For a linear projective space, including a single point, its dual is empty.
All three functions then return `nothing`, as they do for an empty zero locus.
A reduced zero-dimensional locus of two or more points has a hypersurface dual
whose degree is its number of points.

## The computation

Let ``X_k`` be a general linear section of codimension ``k`` in an ``n``-dimensional
polarized variety ``X``, and let ``e_k=\chi_{\mathrm{top}}(X_k)``.
Put ``e_{n+1}=e_{n+2}=0``.
The class formula applied to these sections gives the polar degrees

```math
\begin{equation}
a_k=(-1)^{n-k}(e_k-2e_{k+1}+e_{k+2})
    =\int_X c_{n-k}(J_1(L))c_1(L)^k.
\end{equation}
```

For a nonempty dual, the first nonzero index is its defect, the codimension
is one larger, and the corresponding polar degree is the degree of the dual.
A linear space has only the final value ``a_n=1``; this is handled as an empty
dual instead of a degree-one dual.

These are **topological** Euler characteristics, obtained from
`euler_characteristic(X_k)`, rather than the holomorphic Euler characteristic
of a structure sheaf.
For a zero locus ``Z(E)``, the successive sections are
`zero_locus(E + k * L)`, where `k * L` is the sum of `k` copies of `L`.
For example, the first two are `zero_locus(E + L)` and `zero_locus(E + L + L)`.
The computation stops as soon as the dual invariants are known.

See Tevelev's [*Projective Duality and Homogeneous Spaces*](https://doi.org/10.1007/b138367),
Theorem 7.1 and Example 7.5 in the [arXiv version](https://arxiv.org/abs/math/0112028).
The references below use arXiv numbering.

## Examples from Tevelev

The regression tests cover more than sixty embeddings from these families,
including both hypersurface and defective duals:

| Variety and polarization | Dual degree | Defect | Reference |
|:--|--:|--:|:--|
| Veronese ``\mathbb P^n``, ``\mathcal O(d)``, ``d\geq2`` | ``(n+1)(d-1)^n`` | 0 | Example 7.9 |
| Plane curve of degree ``d``, ``\mathcal O(r)`` | ``d(d+2r-3)`` | 0 | Example 7.6 |
| Smooth curve of genus ``g`` and degree ``d`` | ``2g-2+2d`` | 0 | Example 7.7 |
| Smooth hypersurface of degree ``d\geq2`` in ``\mathbb P^{n+1}`` | ``d(d-1)^n`` | 0 | Example 7.4 |
| Quartic K3 surface, ``\mathcal O(2)`` | 72 | 0 | Example 7.8 |
| Segre ``\mathbb P^m\times\mathbb P^n``, ``1\leq m\leq n`` | ``\binom{n+1}{m}`` | ``n-m`` | Example 2.18 |
| Segre ``\mathbb P^1\times\mathbb P^1\times\mathbb P^1`` | 4 | 0 | §6.2, §7.2 |
| Smooth quadric | 2 | 0 | Example 2.19 |
| ``\mathrm{Gr}(2,2m)``, Plücker embedding | ``m`` | 0 | Example 2.20 |
| ``\mathrm{Gr}(2,5)``, Plücker embedding | 5 | 2 | Example 2.20 |
| Spinor tenfold ``\mathrm{OGr}_+(5,10)`` | 12 | 4 | Example 2.21 |
| Cayley plane | 3 | 0 | Example 2.22 |

`quadric(1)` is represented as a projective line with minimal polarization;
use `line_bundle(quadric(1), 2)` for its conic embedding.

Example 7.18 also gives a formula for the full flag variety of ``\mathrm{SL}_3``
with polarization of Picard degrees ``(a,b)``:

```math
\begin{equation}
\deg X^\vee=12(a^2b+ab^2)-6(a^2+4ab+b^2)+12(a+b)-6.
\end{equation}
```

The tests check several pairs of degrees, as well as the drop in defect under
general hyperplane sections from Theorem 1.30.
They also reproduce the four degrees in issue [#38](https://github.com/HomogeneousTools/PartialFlagVarieties.jl/issues/38):
120 for `Gr(3, 9)`, 126 for `Gr(4, 8)`, 240 for `OGr(8, 16)`, and 72 for `LGr(4, 8)`.

## API

```@docs
projective_dual_degree
projective_dual_codimension
projective_dual_defect
```
