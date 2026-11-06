# Exceptional Collections

Exceptional objects, exceptional sequences, and standard exceptional collections
on partial flag varieties and their zero loci.

## Exceptionality predicates

An object ``\mathcal{E}`` on ``\mathrm{G}/\mathrm{P}`` is **exceptional** if ``\operatorname{RHom}(\mathcal{E},\mathcal{E}) \cong \mathbb{C}``
(concentrated in degree 0). A sequence ``\langle \mathcal{E}_1, \ldots, \mathcal{E}_r \rangle`` is an
**exceptional sequence** if each ``\mathcal{E}_i`` is exceptional and
``\operatorname{RHom}(\mathcal{E}_j, \mathcal{E}_i) = 0`` for all ``i < j``.

```@docs
is_exceptional
is_exceptional_pair
is_strong_exceptional_pair
is_exceptional_sequence
is_strong_exceptional_sequence
is_full_exceptional_sequence
```

## Bundles on zero loci

Define bundles on a zero locus once and use the same exceptionality predicates.
Their base is inferred with [`variety`](@ref); there is no need to pass the
zero locus as an additional argument.

```jldoctest
julia> using PartialFlagVarieties

julia> G = Gr(3, 7);

julia> U = universal_subbundle(G);

julia> A = exterior_power(dual(U), 2);

julia> X = zero_locus(A ⊕ A);

julia> UX = restrict(X, U); UdX = dual(UX); OX = structure_sheaf(X);

julia> is_exceptional(UdX)
true

julia> is_exceptional_pair(OX, UdX)
true

julia> is_strong_exceptional_sequence([UX, OX, UdX])
true
```

The objects can be restrictions, intrinsic bundles such as
[`tangent_bundle(X)`](@ref tangent_bundle(::ZeroLocus)), or their tensor, dual,
exterior and symmetric constructions.
Every pair or sequence must live on the same recorded zero locus.
The sequence methods also accept views and other abstract vectors of
`ZeroLocusBundle`s.

!!! note "Symbolic cohomology"
    These predicates return `true` only when the required Ext dimensions are
    certified. Undetermined cohomology emits a warning, and a `false` result
    in that case does not disprove mathematical exceptionality.
    For example, the current constraints leave one parameter in the self-Ext
    groups of the tangent bundle of a hyperplane ``\mathbb{P}^2 \subset \mathbb{P}^3``,
    although that bundle is exceptional. The predicate retains this ambiguity.
    Inspect `cohomology(dual(F) ⊗ F)` to see the unresolved dimensions.

Fullness still requires a separate mathematical argument; none of these
predicates proves generation of the derived category.

## Restrictions of ambient bundles

Given bundles ``\mathcal{E}_i`` on the ambient variety ``X = \mathrm{G}/\mathrm{P}`` and a zero locus ``Z \subset X``,
the following methods check exceptionality of the **restrictions** ``\mathcal{E}_i|_Z``, computing
``\operatorname{Ext}^\bullet(\mathcal{E}_i|_Z, \mathcal{E}_j|_Z) = \mathrm{H}^\bullet(Z,\, \mathcal{E}_i^\vee \otimes \mathcal{E}_j|_Z)`` via the
Koszul resolution.
These existing methods restrict the arguments and delegate to the same checks
used for `ZeroLocusBundle`s.

```@docs
is_exceptional(::CompletelyReducibleBundle, ::ZeroLocus)
is_exceptional_pair(::CompletelyReducibleBundle, ::CompletelyReducibleBundle, ::ZeroLocus)
is_strong_exceptional_pair(::CompletelyReducibleBundle, ::CompletelyReducibleBundle, ::ZeroLocus)
is_exceptional_sequence(::Vector{<:CompletelyReducibleBundle}, ::ZeroLocus)
is_strong_exceptional_sequence(::Vector{<:CompletelyReducibleBundle}, ::ZeroLocus)
```

## Beilinson collection on ``\mathbb{P}^n``

Beilinson's theorem: the line bundles ``\langle \mathcal{O}, \mathcal{O}(1), \ldots, \mathcal{O}(n) \rangle``
form a full strong exceptional collection on ``\mathbb{P}^n``.

```@docs
beilinson_collection
beilinson_collection_dual
```

## Kapranov collection on quadrics

Kapranov's theorem: on the ``n``-dimensional quadric ``Q^n``, the spinor bundle(s)
together with ``\mathcal{O}, \mathcal{O}(1), \ldots, \mathcal{O}(n-1)`` form a
full exceptional collection.

```@docs
kapranov_collection
```

## Kapranov–Orlov collection on Grassmannians

On a Grassmannian ``\mathrm{Gr}(k,n)``, the Schur functors ``\Sigma^\alpha \mathcal{U}^\vee``
for all partitions ``\alpha`` fitting in a ``k \times (n-k)`` box form a full strong
exceptional collection (Kapranov, 1988).

```@docs
schur_functor
kapranov_bundles_grassmannian
```
