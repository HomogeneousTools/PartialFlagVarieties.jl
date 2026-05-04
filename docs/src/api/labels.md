# Labels

Compact serialization of ambient varieties, bundles, and zero loci via the
`ZeroLocus62` label format.

This interface is mainly useful when you need to exchange data with external
classifications or datasets. It is **not** the primary way to construct
varieties in interactive use; for that, prefer the named constructors or
[`partial_flag_variety`](@ref).

## What is encoded

- `zerolocus62_label(X)` records the ambient homogeneous space.
- `zerolocus62_label(E)` records the ambient space together with the bundle
  summands.
- `zerolocus62_label(Z)` records the ambient space and the defining bundle of
  the zero locus.

Encoding a bundle or zero locus is purely a matter of serialization. It does
**not** certify that the bundle admits a regular section.

## Decoding caveat

The two string-based constructors in the package do different things:

- `partial_flag_variety("A3", [2])` parses a Dynkin type string and marked nodes.
- `PartialFlagVariety("31")` decodes a ZeroLocus62 label.

The one-argument `PartialFlagVariety(::AbstractString)` constructor therefore
belongs conceptually to the label interface, not to the ordinary constructor
family.

## API

```@docs
zerolocus62_label
PartialFlagVariety(::AbstractString)
```
