# Conventions & Notation

This page records the package conventions that most often cause confusion.

## Marked nodes

The package uses the convention

> **marked nodes = crossed-out nodes = nonparabolic nodes**

So `marked_nodes(X)` lists the simple roots omitted from the Levi diagram.
Equivalently, for ``X = \mathrm{G}/\mathrm{P}_I``, the marked set ``I`` is the set of crossed-out
nodes in the Dynkin diagram.

This is the opposite of the convention used in some parts of the literature, so
it is worth keeping in mind whenever you translate formulas or diagrams.

## Dynkin numbering

Dynkin types use the Bourbaki / Oscar numbering conventions. In particular:

- roots are numbered in Bourbaki order,
- type ``\mathrm{E}`` uses the standard Bourbaki branching convention,
- string constructors such as `partial_flag_variety("E6", 1)` follow that numbering.

## Roots, weights, and Picard coordinates

- **Roots** are stored in the simple-root basis.
- **Weights** are stored in the fundamental-weight basis.
- The Picard basis is ordered by the marked nodes returned by [`marked_nodes`](@ref).

This matters whenever you interpret vectors of coefficients:

```julia
X = partial_flag_variety("A3", [1, 3])
marked_nodes(X)          # (1, 3)
anticanonical_degrees(X) # coefficients in the basis [ω₁, ω₃]
```

## Semisimplified versus filtered bundles

The main bundle type is [`CompletelyReducibleBundle`](@ref). It stores the
**semisimplification** of an equivariant bundle: enough for tensor algebra,
Borel–Weil–Bott, and most geometric invariants computed here.

When the filtration itself matters, use [`FilteredBundle`](@ref). This is
especially relevant for the filtered tangent bundle and for the conormal
filtration used in Hodge computations.

## Bundle shorthands

The most common bundles have single-letter shorthands, used throughout the
examples and these guides:

| Shorthand | Long form | Bundle |
|:----------|:----------|:-------|
| `O(X)` | [`structure_sheaf`](@ref) | structure sheaf ``\mathcal{O}_X`` |
| `O(X, d)`, `O(X, [d₁, …])` | [`line_bundle`](@ref) | line bundle ``\mathcal{O}(d)`` |
| `S(X)` | [`universal_subbundle`](@ref) | universal subbundle ``\mathcal{S}`` |
| `Q(X)` | [`universal_quotient_bundle`](@ref) | universal quotient ``\mathcal{Q}`` |
| `T(X)` | [`tangent_bundle`](@ref) | tangent bundle ``\mathrm{T}_X`` |
| `E(X, weights)` | [`CompletelyReducibleBundle`](@ref) | bundle from highest weights |

`S`, `Q`, `T`, and `E` are exported functions, so they occupy those names after
`using PartialFlagVarieties`. You can still reuse the letter as a local variable
inside a function or `let` block — a local binding shadows the function there.

## Indexing conventions

Two indexing conventions show up repeatedly:

1. [`Cohomology`](@ref) objects are **0-indexed**: `H[0]` means ``\mathrm{H}^0``.
2. Hodge matrices use ordinary Julia matrix indexing, so the entry
   corresponding to ``(p,q)`` is stored at `[p+1, q+1]`.
3. The Hochschild [`PolyvectorParallelogram`](@ref) is **0-indexed**:
   `P[p, q]` is ``\mathrm{h}^q(X, \bigwedge\nolimits^p \mathrm{T}_X)`` and `P[n]` is
   ``\dim \mathrm{HH}^n``.

For example:

```julia
X = projective_space(2)
H = hodge_numbers(X)
H[2, 2]  # h^{1,1}
```

## Zero-locus assumptions

`zero_locus(E)` studies the zero locus of a **regular section** of `E`. The
package checks basic consistency conditions such as `rank(E) <= dimension(X)`,
but it does **not** prove that a regular section exists or that the geometric
zero locus is nonempty.

This assumption is used throughout the zero-locus, Hodge, and Hochschild APIs.

## Integer and degree-vector twists

Functions such as [`twisted_hodge_numbers`](@ref) accept an integer twist
when the Picard group is rank 1 (`j` means tensoring by the generator
normalized by the marked node) and a **degree vector** with one entry per
marked node otherwise, mirroring the `O(X, [d₁, …])` constructor:

```julia
X = projective_space(4)
twisted_hodge_numbers(X, 1)        # Ω^p_X(1)

Y = partial_flag_variety("A3", [1, 3])
twisted_hodge_numbers(Y, [1, 2])   # Ω^p_Y ⊗ O(1, 2)
```

An explicit line bundle works as well in either case.

## Dynkin strings versus ZeroLocus62 labels

Two string syntaxes coexist:

- `partial_flag_variety("A3", [2])` parses a **Dynkin type string** plus marked nodes.
- `PartialFlagVariety("31")` decodes a **ZeroLocus62 label**.

They serve different purposes. The first is for ordinary user input; the second
is for compact serialization and dataset interoperability; see [Labels](api/labels.md).
