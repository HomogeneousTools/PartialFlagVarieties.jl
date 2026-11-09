# Numerical invariants

`volume(X)` computes the anticanonical top self-intersection from the exact
symbolic Hilbert polynomial. For a zero locus it requires the anticanonical
class to lift to an ample line bundle on the ambient flag variety. Failure
of this sufficient condition leaves its Fano status undecided.

`invariants(X)` collects dimension, anticanonical volume, holomorphic Euler
characteristics, and Hodge numbers in a named tuple. Set `hodge=false` to
skip Hodge computation, and `max_tangent_power=k` to choose how many exterior
tangent powers to include. A volume not certified by the ambient criterion
is returned as `nothing`.

```@docs
volume
anticanonical_volume
invariants
```
