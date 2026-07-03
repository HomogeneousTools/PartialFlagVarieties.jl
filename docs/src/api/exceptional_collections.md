# Exceptional Collections

Exceptional objects, exceptional sequences, and standard exceptional collections
on partial flag varieties and their zero loci.

## Exceptionality predicates

An object ``\mathcal{E}`` on ``\mathrm{G}/\mathrm{P}`` is **exceptional** if ``\operatorname{RHom}(\mathcal{E},\mathcal{E}) \cong \mathbb{C}``
(concentrated in degree 0). A sequence ``\langle \mathrm{E}_1, \ldots, \mathcal{E}_r \rangle`` is an
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

## Exceptionality on zero loci

Given bundles ``\mathcal{E}_i`` on the ambient variety ``X = \mathrm{G}/\mathrm{P}`` and a zero locus ``Z \subset X``,
the following methods check exceptionality of the **restrictions** ``\mathcal{E}_i|_Z``, computing
``\operatorname{Ext}^\vee(\mathcal{E}_i|_Z, \mathcal{E}_j|_Z) = \mathrm{H}^\bullet(Z,\, \mathcal{E}_i^\vee \otimes \mathcal{E}_j|_Z)`` via the
Koszul resolution.

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
