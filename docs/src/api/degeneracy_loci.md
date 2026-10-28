# Degeneracy loci

For a map $E\to F$ of ranks $e,f$, the locus of rank at most~$r$ has expected
codimension~$(e-r)(f-r)$.
On~$\operatorname{Gr}_X(e-r,E)$ it has the kernel incidence presentation
as a zero locus of~$S^\vee\otimes\pi^*F$.
The projection is an isomorphism when the lower rank stratum is absent.
`check_smoothness=true` verifies global generation and requires the lower
stratum's expected codimension to exceed the base dimension.
`degeneracy_locus` returns a `DegeneracyLocus` retaining the map data and its
incidence resolution. `determinantal_resolution` returns the resolution directly.
Without the smoothness check, the image can be singular and its invariants
can differ from those computed on the resolution.
Regularity and expected dimension remain assumptions about the general section.

```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(8);

julia> Z = degeneracy_locus(4 * structure_sheaf(X), 4 * line_bundle(X, 1), 2);

julia> dimension(Z)
4

julia> Z isa DegeneracyLocus
true

julia> incidence_resolution(Z) isa ZeroLocus  # presentation on P⁸ × Gr(2,4)
true

julia> R = degeneracy_locus(4 * structure_sheaf(X), 4 * line_bundle(X, 1), 2;
                          presentation=:relative);

julia> cohomology(structure_sheaf(R))[0]
1
```

When~$E=L^{\oplus e}$, the Grassmann bundle is a product and its defining bundle
is~$S_0^\vee\boxtimes(F\otimes L^\vee)$, up to the order of the factors.
The default constructor recognizes this case, including after dualizing the map.
`presentation=:relative` forces the general algorithm.
Supply `base_equations=B` when the base is itself the regular zero locus of~$B$
in the homogeneous variety carrying `E` and `F`.

On a degeneracy locus, use `restrict`, `cotangent_bundle`, `tangent_bundle`,
`exterior_power`, `canonical_bundle`, `cohomology`, `hodge_numbers`, and
`hochschild_cohomology`.
`chi(bundle)` is exact even if cohomology is symbolic.
`degree(anticanonical_bundle(Z))` computes the anticanonical top intersection.
The Koszul resolution, fibrewise Borel--Weil--Bott, and the conormal sequence
reduce these calculations to cohomology on the homogeneous base.
Unknown ranks of induced maps remain symbolic; the code does not assume
maximal rank or degeneration of an arbitrary spectral sequence.
`ambient_variety(Z)` is the homogeneous base containing the degeneracy image.
Use `ambient_variety(incidence_resolution(Z))` to inspect the auxiliary ambient.
Bundles on the locus retain it as their `variety`; their calculations use the
underlying incidence presentation.

## Pfaffian and orbital constructions

`pfaffian_locus(E,L; corank=c)` handles sections of~$\bigwedge^2E\otimes L$
of corank at least~$c$, with expected codimension~$c(c-1)/2$.
The allowed rank~$\operatorname{rank}(E)-c$ must be even; the default is~$c=3$.
The incidence equations form the quotient
$\bigwedge^2\pi^*E/\bigwedge^2S$, tensored by~$\pi^*L$.
Its filtration is retained, since replacing it by a direct sum would give
incorrect cohomology in general.

`symmetric_locus(E,L; corank=c)` uses the filtered quotient
$\operatorname{Sym}^2\pi^*E/\operatorname{Sym}^2S$, tensored by~$\pi^*L$.
Its expected codimension is~$c(c+1)/2$, and the default is~$c=1$.
Both constructors support even and odd bundle ranks, arbitrary admissible
coranks, and filtered input and base equations.

## Optional geometry checks

The determinantal, Pfaffian, symmetric and orbital locus constructors accept
`check_global_generation=false`, `check_smoothness=false`, and `check_nonempty=false`.
Global generation is exact for completely reducible homogeneous section bundles.
For filtered bundles, dominant graded pieces give a sufficient certificate;
failure of that certificate does not prove the extension is not generated.
The smoothness check also verifies global generation of the base equations
and excludes the deeper stratum by expected codimension.
The nonemptiness check computes the top intersection of a certified ample
polarization on the expected regular locus; it can be expensive.
A failed check raises `ArgumentError`.
These checks do not choose a section or prove regularity for a particular section.
`determinantal_resolution` accepts generation and nonemptiness checks separately.

## Anticanonical positivity

Use `is_strongly_fano(Z)` to query the sufficient ambient ampleness certificate
after constructing a locus.
The constructors return `DegeneracyLocus`, whose incidence presentation is
a `ZeroLocus` or `RelativeZeroLocus`.
The query uses the saved deeper-stratum bound to identify the image with
its resolution, then looks for an ample anticanonical lift on the incidence
ambient or an ample anticanonical line descending from the base.
This can certify a Fano degeneracy locus even when
`is_strongly_fano(incidence_resolution(Z))` is false.
A negative result does not prove that the locus is not Fano:
the restriction can be ample even when the ambient lift is not.
These queries assume a smooth, nonempty regular locus and do not certify
those assumptions themselves.
No general `is_fano` predicate is provided.

`partially_decomposable_locus(E)` implements the rank-six three-form construction
of [Benedetti–Filippini–Manivel–Tanturri](https://arxiv.org/abs/1704.01436).
Its equations are~$\bigwedge^3Q$ on~$\mathbb P(E)$.
Appendix A is therefore directly within this API.
`rank_one_nilpotent_locus(E,L)` implements the rank-one traceless endomorphism
construction; this is a specific orbital family, not a general orbit-closure engine.

## Resolutions on the base

`eagon_northcott(E,F)` resolves the structure sheaf of the **maximal-minor** locus.
`buchsbaum_rim(E,F)` resolves the **cokernel** of the map, a different sheaf.
Both require~$\operatorname{rank}E\geq\operatorname{rank}F>0$ and the expected grade.
Neither is used as a resolution for arbitrary smaller minors.
`gulliksen_negard(E,F)` resolves the **submaximal-minor** locus of a square map
under the expected codimension-four hypothesis.
`pfaffian_resolution(E,L)` is the codimension-three Buchsbaum–Eisenbud resolution.
Twist a resolution with `tensor_product(R,L)` and then use `chi` or `cohomology`.
These direct resolutions can be much smaller than the incidence Koszul complex.

## Computational limits and examples

The incidence method supports filtered homogeneous input bundles
and bases cut out by filtered homogeneous bundles.
High ranks can make Schur decompositions and Koszul powers expensive, and
individual cohomology groups can remain symbolic even when Euler invariants are exact.

All five examples in Table 8 of
[arXiv:1704.01436](https://arxiv.org/abs/1704.01436) are implemented in
`examples/DegeneracyLoci.jl`:

| Row | Anticanonical degree | χ(Ω¹) | χ(Ω²) | h⁰(ω∨) |
| --- | ---: | ---: | ---: | ---: |
| 1 | 63 | -2 | 21 | 19 |
| 2 | 69 | -4 | 26 | 20 |
| 3 | 47 | -7 | 54 | 16 |
| 4 | 40 | -18 | 114 | 15 |
| 5 | 70 | -6 | 46 | 21 |

The first three use the relative projective construction from Appendix A.
The last two have smaller presentations on `Fl(1,3;4)` and `Fl(1,4;5)`:
a quadric section and three hyperplane sections, respectively.
`rank_one_nilpotent_locus` also provides their general incidence construction.

Run `julia --project=. examples/DegeneracyLoci.jl` for these Euler invariants.
After including the example file, `hodge_numbers(orbital_fano_fourfolds()[i])`
computes the possibly symbolic Hodge diamond of row `i`.

```@docs
degeneracy_locus
DegeneracyLocus
DegeneracyLocusBundle
incidence_resolution
determinantal_resolution
RelativeZeroLocus
RelativeZeroLocusBundle
pfaffian_locus
symmetric_locus
partially_decomposable_locus
rank_one_nilpotent_locus
BundleResolution
eagon_northcott
gulliksen_negard
buchsbaum_rim
pfaffian_resolution
```
