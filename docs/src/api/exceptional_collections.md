# Exceptional Collections

Exceptional objects, exceptional sequences, and standard exceptional collections
on partial flag varieties.

## Exceptionality predicates

An object ``E`` on ``G/P`` is **exceptional** if ``\operatorname{RHom}(E,E) \cong \mathbb{k}``
(concentrated in degree 0). A sequence ``\langle E_1, \ldots, E_r \rangle`` is an
**exceptional sequence** if each ``E_i`` is exceptional and
``\operatorname{RHom}(E_j, E_i) = 0`` for all ``i < j``.

```@docs
is_exceptional
is_exceptional_pair
is_strong_exceptional_pair
is_exceptional_sequence
is_strong_exceptional_sequence
is_full_exceptional_sequence
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
