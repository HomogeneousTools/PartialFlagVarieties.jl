# ═══════════════════════════════════════════════════════════════════════════════
#  Constructions — named convenience constructors for classical varieties
#
#  Provides familiar notation for well-known partial flag varieties:
#  Grassmannians, projective spaces, quadrics, isotropic Grassmannians,
#  exceptional varieties, etc.
# ═══════════════════════════════════════════════════════════════════════════════

export Gr, OGr, SGr, LGr, IGr
export projective_space, quadric
export cayley_plane, freudenthal_variety
export adjoint_variety, coadjoint_variety
export flag_variety

# ═══════════════════════════════════════════════════════════════════════════════
#  Type A: Grassmannians and flags
# ═══════════════════════════════════════════════════════════════════════════════

"""
    Gr(k, n) -> PartialFlagVariety

The Grassmannian ``\\mathrm{Gr}(k, n)`` of ``k``-planes in ``\\mathbb{C}^n``.

It is the generalized Grassmannian ``A_{n-1}/P_k``, encoded by marking node
``k`` in the Dynkin diagram of type ``A_{n-1}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = Gr(2, 5);

julia> dimension(V)
6

julia> euler_characteristic(V)
10
```
"""
function Gr(k::Integer, n::Integer)
  k, n = Int(k), Int(n)
  1 <= k <= n - 1 || throw(ArgumentError("Gr($k, $n): need 1 ≤ k ≤ n-1"))
  DT = TypeA{n - 1}
  return partial_flag_variety(DT, (k,), "Gr($k, $n)")
end

"""
    projective_space(n) -> PartialFlagVariety

The projective space ``\\mathbb{P}^n = \\mathrm{Gr}(1, n+1)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = projective_space(4);

julia> dimension(V)
4

julia> betti_numbers(V)
5-element Vector{BigInt}:
 1
 1
 1
 1
 1
```
"""
function projective_space(n::Integer)
  n = Int(n)
  n >= 1 || throw(ArgumentError("projective_space($n): need n ≥ 1"))
  DT = TypeA{n}
  return partial_flag_variety(DT, (1,), "ℙ$n")
end

"""
    flag_variety(n, dimensions) -> PartialFlagVariety

The type-``A`` partial flag variety
``\\mathrm{Fl}(d_1, \\ldots, d_r; n) = A_{n-1}/P_{\\{d_1,\\ldots,d_r\\}}``.

Here `dimensions = (d_1, ..., d_r)` are the dimensions of the subspaces in the
partial flag
``0 \\subsetneq V_{d_1} \\subsetneq \\cdots \\subsetneq V_{d_r} \\subsetneq \\mathbb{C}^n``;
they must be distinct, sorted in increasing order, and satisfy ``1 \\leq d_i \\leq n-1``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = flag_variety(4, [1, 2]);

julia> dimension(V)
5
```
"""
function flag_variety(n::Integer, dimensions::Vector{<:Integer})
  n = Int(n)
  all(1 <= d <= n - 1 for d in dimensions) ||
    throw(ArgumentError("flag_variety($n, $dimensions): need 1 ≤ dᵢ ≤ $(n - 1)"))
  issorted(dimensions) ||
    throw(ArgumentError("flag_variety($n, $dimensions): dimensions must be sorted"))
  allunique(dimensions) ||
    throw(ArgumentError("flag_variety($n, $dimensions): dimensions must be distinct"))
  flag = join(dimensions, ",")
  return partial_flag_variety(TypeA{n - 1}, Tuple(Int.(dimensions)), "Fl($(flag); $n)")
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Type B/D: Orthogonal Grassmannians
# ═══════════════════════════════════════════════════════════════════════════════

"""
    OGr(k, n) -> PartialFlagVariety

The orthogonal Grassmannian ``\\mathrm{OGr}(k, n)`` of isotropic ``k``-planes
in ``\\mathbb{C}^n`` with a symmetric bilinear form.

For ``n = 2m + 1`` (odd): type ``B_m``, mark node ``k``.
For ``n = 2m`` (even) and ``k < m``: type ``D_m``, mark node ``k``.
For ``n = 2m`` and ``k = m``: this is the spinor variety ``\\mathrm{OGr}_+(m, 2m)``,
type ``D_m`` mark node ``m``.

When ``n = 2m`` and ``k = m``, the maximal orthogonal Grassmannian has two
connected components. `OGr(m, 2m)` picks the component corresponding to node
``m``; the other component is ``D_m/P_{m-1}``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = OGr(5, 10);  # spinor variety S₅

julia> dimension(V)
10

julia> euler_characteristic(V)
16
```
"""
function OGr(k::Integer, n::Integer)
  k, n = Int(k), Int(n)
  1 <= k || throw(ArgumentError("OGr($k, $n): need k ≥ 1"))
  if isodd(n)
    m = (n - 1) ÷ 2
    1 <= k <= m || throw(ArgumentError("OGr($k, $n): need k ≤ $m for B$m"))
    DT = TypeB{m}
    return partial_flag_variety(DT, (k,), "OGr($k, $n)")
  else
    m = n ÷ 2
    1 <= k <= m || throw(ArgumentError("OGr($k, $n): need k ≤ $m for D$m"))
    DT = TypeD{m}
    return partial_flag_variety(DT, (k,), "OGr($k, $n)")
  end
end

"""
    SGr(k, n) -> PartialFlagVariety

The symplectic Grassmannian ``\\mathrm{SGr}(k, n)`` of isotropic ``k``-planes
in ``\\mathbb{C}^n`` with a skew-symmetric form (``n`` must be even).

Type ``C_{n/2}``, mark node ``k``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = SGr(2, 6);

julia> dimension(V)
7
```
"""
function SGr(k::Integer, n::Integer)
  k, n = Int(k), Int(n)
  iseven(n) || throw(ArgumentError("SGr($k, $n): n must be even"))
  m = n ÷ 2
  1 <= k <= m || throw(ArgumentError("SGr($k, $n): need 1 ≤ k ≤ $m"))
  DT = TypeC{m}
  return partial_flag_variety(DT, (k,), "SGr($k, $n)")
end

"""
    LGr(n, 2n) -> PartialFlagVariety

The Lagrangian Grassmannian ``\\mathrm{LGr}(n, 2n) = \\mathrm{SGr}(n, 2n)`` of
Lagrangian (maximal isotropic) ``n``-planes in a ``2n``-dimensional symplectic
space. The second argument must equal ``2n``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = LGr(3, 6);

julia> dimension(V)
6
```
"""
function LGr(n::Integer, m::Integer)
  n, m = Int(n), Int(m)
  n >= 1 || throw(ArgumentError("LGr($n, $m): need n ≥ 1"))
  m == 2n || throw(ArgumentError("LGr($n, $m): second argument must be 2n = $(2n)"))
  return SGr(n, m)
end

"""
    IGr(k, n) -> PartialFlagVariety

Synonym for [`OGr`](@ref), used to emphasize that the planes are isotropic for
an orthogonal form.
"""
IGr(k::Integer, n::Integer) = OGr(k, n)

"""
    quadric(n) -> PartialFlagVariety

The smooth quadric hypersurface ``Q_n \\subset \\mathbb{P}^{n+1}``.

For odd ``n = 2m - 1``: ``\\mathrm{OGr}(1, 2m + 1)`` = ``B_m / P_1``.
For even ``n = 2m - 2``: ``\\mathrm{OGr}(1, 2m)`` = ``D_m / P_1``.

Equivalently, `quadric(n)` is the isotropic-line Grassmannian
``\\mathrm{OGr}(1, n+2)``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = quadric(4);

julia> dimension(V)
4

julia> betti_numbers(V)
5-element Vector{BigInt}:
 1
 1
 2
 1
 1
```
"""
function quadric(n::Integer)
  n = Int(n)
  n >= 1 || throw(ArgumentError("quadric($n): need n ≥ 1"))
  n == 1 && return partial_flag_variety(TypeA{1}, (1,), "Q$n")
  n == 2 &&
    return partial_flag_variety(ProductDynkinType{Tuple{TypeA{1},TypeA{1}}}, (1, 2), "Q$n")
  if isodd(n)
    m = (n + 1) ÷ 2
    return partial_flag_variety(TypeB{m}, (1,), "Q$n")
  else
    m = n ÷ 2 + 1
    return partial_flag_variety(TypeD{m}, (1,), "Q$n")
  end
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Exceptional varieties
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cayley_plane() -> PartialFlagVariety

The Cayley plane ``\\mathbb{OP}^2 = E_6 / P_1``, the 16-dimensional
cominuscule variety.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = cayley_plane();

julia> dimension(V)
16

julia> euler_characteristic(V)
27
```
"""
function cayley_plane()
  return partial_flag_variety(TypeE{6}, (1,), "OP²")
end

"""
    freudenthal_variety() -> PartialFlagVariety

The Freudenthal variety ``E_7 / P_7``, the 27-dimensional cominuscule variety.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = freudenthal_variety();

julia> dimension(V)
27

julia> euler_characteristic(V)
56
```
"""
function freudenthal_variety()
  return partial_flag_variety(TypeE{7}, (7,), "E₇/P₇")
end

"""
    adjoint_variety(::Type{DT}) -> PartialFlagVariety

The adjoint variety of a **simple** Dynkin type `DT`, i.e. the projectivisation
of the minimal nilpotent orbit.

The marked node is the adjoint node, determined by the highest root of `DT`.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = adjoint_variety(TypeE{6});

julia> dimension(V)
21
```
"""
function adjoint_variety(::Type{DT}) where {DT<:SimpleDynkinType}
  # Determine the adjoint node
  DT <: TypeA && return _adjoint_A(DT)
  DT <: TypeB && return partial_flag_variety(DT, (2,), "Adj($(Semisimple._type_name(DT)))")
  DT <: TypeC && return partial_flag_variety(DT, (1,), "Adj($(Semisimple._type_name(DT)))")
  DT <: TypeD && return partial_flag_variety(DT, (2,), "Adj($(Semisimple._type_name(DT)))")
  DT <: TypeE{6} && return partial_flag_variety(DT, (2,), "Adj(E₆)")
  DT <: TypeE{7} && return partial_flag_variety(DT, (1,), "Adj(E₇)")
  DT <: TypeE{8} && return partial_flag_variety(DT, (8,), "Adj(E₈)")
  DT <: TypeF4 && return partial_flag_variety(DT, (1,), "Adj(F₄)")
  DT <: TypeG2 && return partial_flag_variety(DT, (2,), "Adj(G₂)")

  error("Adjoint variety not implemented for $DT")
end

function _adjoint_A(::Type{DT}) where {DT<:TypeA}
  R = rank(DT)
  return partial_flag_variety(DT, (1, R), "Adj($(Semisimple._type_name(DT)))")
end

"""
    coadjoint_variety(::Type{DT}) -> PartialFlagVariety

The coadjoint variety of a **simple** Dynkin type `DT`.

For simply laced types it agrees with the adjoint variety. For nonsimply laced
types it is determined by the highest short root.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> V = coadjoint_variety(TypeG2);

julia> dimension(V)
5
```
"""
function coadjoint_variety(::Type{DT}) where {DT<:SimpleDynkinType}
  DT <: TypeA && return _adjoint_A(DT)  # self-dual
  DT <: TypeB &&
    return partial_flag_variety(DT, (1,), "Coadj($(Semisimple._type_name(DT)))")
  DT <: TypeC &&
    return partial_flag_variety(DT, (2,), "Coadj($(Semisimple._type_name(DT)))")
  DT <: TypeD &&
    return partial_flag_variety(DT, (2,), "Coadj($(Semisimple._type_name(DT)))")
  DT <: TypeE{6} && return partial_flag_variety(DT, (2,), "Coadj(E₆)")
  DT <: TypeE{7} && return partial_flag_variety(DT, (1,), "Coadj(E₇)")
  DT <: TypeE{8} && return partial_flag_variety(DT, (8,), "Coadj(E₈)")
  DT <: TypeF4 && return partial_flag_variety(DT, (4,), "Coadj(F₄)")
  DT <: TypeG2 && return partial_flag_variety(DT, (1,), "Coadj(G₂)")

  error("Coadjoint variety not implemented for $DT")
end
