# ═══════════════════════════════════════════════════════════════════════════════
#  FanoThreefolds.jl — Hodge diamonds and polyvector parallelograms
#  of homogeneous Fano threefolds
#
#  Varieties are taken from:
#    P. Belmans, "homogeneous-description.m2"
#    https://github.com/pbelmans/bivector-fields-fano-3-folds
#
#  Each variety is realised as the zero locus Z(s) ⊂ G/P of a regular section
#  of a completely reducible equivariant bundle E on a product of Grassmannians.
#  Labels follow the Mori–Mukai classification (ρ-serial_number).
#
#  For the classical G/P Fano threefolds (ℙ³, Q³, Fl, products) the polyvector
#  parallelogram H^q(X, ∧^p T_X) is also computed.
#
#  Usage:
#    julia --project=. examples/FanoThreefolds.jl
# ═══════════════════════════════════════════════════════════════════════════════

using PartialFlagVarieties
using PrettyTables
using Lie

# ─── GL(n) weight → fundamental weight coordinates ───────────────────────────
#  Identical to FanoFourfolds.jl.

function gl_weight_to_omega_flag(ks::Vector{Int}, n::Int, w::Vector{Int})
  eps = Int[]
  append!(eps, w[n - ks[1] + 1:n])
  for i in 2:length(ks)
    append!(eps, w[n - ks[i] + 1:n - ks[i - 1]])
  end
  append!(eps, w[1:n - ks[end]])
  [eps[i] - eps[i + 1] for i in 1:n - 1]
end

# ─── Build ambient variety from M2 factor list ───────────────────────────────
#  factors = [[k1,n1], [k2,n2], ...] where {ki,ni} = Gr(ki,ni) in M2 notation.

function build_ambient(factors::Vector)
  factor_types = Any[]
  factor_marks = Vector{Int}[]
  for f in factors
    fv = Int.(f)
    n = fv[end]
    ks = fv[1:end-1]
    push!(factor_types, TypeA{n - 1})
    push!(factor_marks, ks)
  end

  if length(factor_types) == 1
    DT = factor_types[1]
    return partial_flag_variety(DT, Tuple(factor_marks[1]))
  end

  DT = factor_types[1]
  offset = rank(DT)
  all_marks = Int[factor_marks[1]...]
  for i in 2:length(factor_types)
    DT = ProductDynkinType{Tuple{DT, factor_types[i]}}
    append!(all_marks, [m + offset for m in factor_marks[i]])
    offset += rank(factor_types[i])
  end
  partial_flag_variety(DT, Tuple(all_marks))
end

# ─── Build equivariant bundle from M2 bundle specification ───────────────────
#  bundle_weights = list of summands; each summand = list of per-factor
#  GL weight vectors (one entry per element of `factors`).

function build_m2_bundle(X, factors, bundle_weights)
  mdt = marked_dynkin_type(X)
  DT = PartialFlagVarieties._ambient_type(mdt)

  factor_ks = [f[1:end-1] for f in factors]
  factor_ns = [f[end] for f in factors]

  summands = IrrepLevi[]
  for summand in bundle_weights
    omega_coords = Int[]
    for (j, w) in enumerate(summand)
      omega = gl_weight_to_omega_flag(factor_ks[j], factor_ns[j], Int.(w))
      append!(omega_coords, omega)
    end
    lam = WeightLatticeElem(DT, omega_coords)
    push!(summands, IrrepLevi(mdt, lam))
  end
  CompletelyReducibleBundle(X, summands)
end

# ─── M2 variety data ─────────────────────────────────────────────────────────
#
#  Format: (MM_label, factors, bundle_summands, description)
#
#  Source: homogeneous-description.m2 by P. Belmans

const M2_FANO3 = [
  # ── ρ = 1, single Grassmannian ─────────────────────────────────────────────
  ("1-5",  [[2,5]],
    [[[0,0,0,2,2]], [[0,0,0,1,1]], [[0,0,0,1,1]]],
    "O(2)⊕O(1)² in Gr(2,5)"),

  ("1-6",  [[2,5]],
    [[[0,0,0,2,1]], [[0,0,0,1,1]]],
    "S_{2,1}U⊕O(1) in Gr(2,5)"),

  ("1-7",  [[2,6]],
    [[[0,0,0,0,1,1]], [[0,0,0,0,1,1]], [[0,0,0,0,1,1]],
     [[0,0,0,0,1,1]], [[0,0,0,0,1,1]]],
    "O(1)⁵ in Gr(2,6)"),

  ("1-8",  [[3,6]],
    [[[0,0,0,1,1,0]], [[0,0,0,1,1,1]], [[0,0,0,1,1,1]], [[0,0,0,1,1,1]]],
    "∧²U ⊕ (det U)³ in Gr(3,6)"),

  ("1-9",  [[2,7]],
    [[[1,0,0,0,0,1,1]], [[0,0,0,0,0,1,1]], [[0,0,0,0,0,1,1]]],
    "(Q⊗det U) ⊕ O(1)² in Gr(2,7)"),

  ("1-10", [[3,7]],
    [[[0,0,0,0,1,1,0]], [[0,0,0,0,1,1,0]], [[0,0,0,0,1,1,0]]],
    "(∧²U)³ in Gr(3,7)"),

  ("1-15", [[2,5]],
    [[[0,0,0,1,1]], [[0,0,0,1,1]], [[0,0,0,1,1]]],
    "O(1)³ in Gr(2,5)"),

  # ── ρ = 2, product of two Grassmannians ────────────────────────────────────
  ("2-14", [[1,2],[2,5]],
    [[[0,1],[0,0,0,1,1]], [[0,0],[0,0,0,1,1]],
     [[0,0],[0,0,0,1,1]], [[0,0],[0,0,0,1,1]]],
    "bundle on ℙ¹×Gr(2,5)"),

  ("2-17", [[2,4],[1,4]],
    [[[0,0,1,0],[0,0,0,1]], [[0,0,1,1],[0,0,0,1]], [[0,0,1,1],[0,0,0,0]]],
    "bundle on Gr(2,4)×ℙ³"),

  ("2-20", [[1,3],[2,5]],
    [[[0,0,1],[0,0,0,1,0]], [[0,0,0],[0,0,0,1,1]],
     [[0,0,0],[0,0,0,1,1]], [[0,0,0],[0,0,0,1,1]]],
    "bundle on ℙ²×Gr(2,5)"),

  ("2-21", [[1,5],[2,4]],
    [[[0,0,0,0,1],[0,0,1,0]], [[0,0,0,0,0],[0,0,1,1]],
     [[0,0,0,0,1],[0,0,1,0]]],
    "bundle on ℙ⁴×Gr(2,4)"),

  ("2-22", [[1,4],[2,5]],
    [[[0,0,0,1],[1,1,0,1,1]], [[0,0,0,0],[0,0,0,1,1]],
     [[0,0,0,0],[0,0,0,1,1]], [[0,0,0,0],[0,0,0,1,1]]],
    "bundle on ℙ³×Gr(2,5)"),

  ("2-26", [[2,4],[2,5]],
    [[[1,0,1,1],[0,0,0,1,0]], [[0,0,1,1],[0,0,0,0,0]],
     [[0,0,0,0],[0,0,0,1,1]], [[0,0,0,0],[0,0,0,1,1]]],
    "bundle on Gr(2,4)×Gr(2,5)"),

  # ── ρ = 9 ──────────────────────────────────────────────────────────────────
  ("9-1",  [[1,2],[1,3],[1,2]],
    [[[0,2],[0,0,2],[0,0]]],
    "O(-2)⊠O(-2)⊠O on ℙ¹×ℙ²×ℙ¹"),
]

# ─── Homogeneous G/P Fano threefolds (for polyvector parallelogram) ───────────

const GP_FANO3 = [
  ("ℙ³",
    () -> projective_space(3),
    "A₃/P₁, index 4"),
  ("Q³",
    () -> quadric(3),
    "B₂/P₁, index 3"),
  ("Fl(1,2;3)",
    () -> partial_flag_variety(TypeA{2}, (1, 2)),
    "A₂/P_{1,2}"),
  ("ℙ¹×ℙ²",
    () -> partial_flag_variety(
      ProductDynkinType{Tuple{TypeA{1},TypeA{2}}}, (1, 2)),
    "A₁×A₂/P_{1,2}"),
  ("ℙ¹×ℙ¹×ℙ¹",
    () -> partial_flag_variety(
      ProductDynkinType{Tuple{TypeA{1},TypeA{1},TypeA{1}}}, (1, 2, 3)),
    "A₁³/P_{1,2,3}"),
]

# ─── Hodge diamond display ───────────────────────────────────────────────────

"""
Print a Hodge diamond for a threefold in the traditional centered layout.

h[p+1, q+1] = h^{p,q}.  Within each anti-diagonal n = p+q, entries are
printed with p decreasing (h^{n,0} on the left, h^{0,n} on the right).
Rows are indented by `ceil((w+1)/2)` characters per step from the widest
row, so the diamond is visually centred.
"""
# ─── Main computation ────────────────────────────────────────────────────────

function main()
  println()
  println("=" ^ 72)
  println("  Homogeneous Fano threefolds (Belmans M2 source)")
  println("=" ^ 72)

  # ── Zero-locus varieties from M2 file ──────────────────────────────────────
  zl_results = []
  zl_errors  = String[]

  for (label, factors, bweights, desc) in M2_FANO3
    local X, E, Z
    try
      X = build_ambient(factors)
      E = build_m2_bundle(X, factors, bweights)
      Z = zero_locus(E)
    catch err
      push!(zl_errors, "$label: construction failed — $err")
      continue
    end
    d = dimension(Z)
    if d != 3
      push!(zl_errors, "$label: expected dim 3, got $d")
      continue
    end
    h = hodge_numbers(Z)
    χ = euler_characteristic(Z)
    push!(zl_results, (; label, Z, d, h, χ, desc))
  end

  for e in zl_errors
    println("ERROR: $e")
  end

  # ── Hodge diamonds (zero-locus varieties) ──────────────────────────────────
  println("\n─── Hodge diamonds (M2 zero-locus varieties) ───\n")
  for r in zl_results
    println("  $(r.label)  $(r.desc)  (χ=$(r.χ))")
    print_hodge_diamond(stdout, r.h)
    println()
  end

  # ── Summary table ──────────────────────────────────────────────────────────
  println("─── Summary ───\n")
  header = ["Label", "Description", "χ", "h¹¹", "h¹²", "h²¹", "h²²"]
  tdata  = Matrix{Any}(undef, length(zl_results), length(header))
  for (i, r) in enumerate(zl_results)
    h = r.h
    tdata[i, 1] = r.label
    tdata[i, 2] = r.desc
    tdata[i, 3] = r.χ
    tdata[i, 4] = h[2, 2]   # h^{1,1}
    tdata[i, 5] = h[2, 3]   # h^{1,2}
    tdata[i, 6] = h[3, 2]   # h^{2,1}
    tdata[i, 7] = h[3, 3]   # h^{2,2}
  end
  pretty_table(tdata;
    column_labels=header,
    alignment=[:l, :l, :c, :c, :c, :c, :c],
  )

  # ── G/P varieties: polyvector parallelograms ───────────────────────────────
  println("\n─── G/P Fano threefolds: polyvector parallelograms ───\n")
  for (name, ctor, desc) in GP_FANO3
    X = ctor()
    h = hodge_numbers(X)
    pvp = hochschild_cohomology(X)
    println("  $name  ($desc,  χ=$(euler_characteristic(X)))")
    print_hodge_diamond(stdout, h)
    println()
    show(stdout, MIME"text/plain"(), pvp)
    println("\n")
  end

  println()
end

main()
