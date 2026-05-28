# ═══════════════════════════════════════════════════════════════════════════════
#  Semisimple.jl — extensions for inclusion in the Semisimple.jl package
#
#  Provides:
#   - cartan_type_with_ordering: identify a Cartan matrix → Dynkin type + ordering
#   - cartan_type: identify a Cartan matrix → Dynkin type
#   - _cartan_type_to_dynkin_type: convert [(Symbol, Int), ...] → DynkinType
#   - DynkinType(::AbstractString): parse "A3", "A2xB3", etc.
#
#  The algorithm for cartan_type_with_ordering follows the approach from
#  OSCAR.jl (Oscar.jl/src/LieTheory/CartanMatrix.jl).
# ═══════════════════════════════════════════════════════════════════════════════

export cartan_type, cartan_type_with_ordering
export parse_dynkin_type

"""
  Semisimple.WeightLatticeElem(::Type{DT}) -> WeightLatticeElem{DT,R}

Construct the zero weight in the weight lattice of `DT`.

# Examples
```julia
julia> WeightLatticeElem(TypeA{2})
0
```
"""
function Semisimple.WeightLatticeElem(::Type{DT}) where {DT<:DynkinType}
  R = rank(DT)
  return WeightLatticeElem(DT, zero(SVector{R,Int}))
end

Base.zero(::Type{Semisimple.WeightLatticeElem{DT,R}}) where {DT<:DynkinType,R} =
  WeightLatticeElem(
    DT
  )

# ═══════════════════════════════════════════════════════════════════════════════
#  cartan_type_with_ordering
# ═══════════════════════════════════════════════════════════════════════════════

"""
    cartan_type_with_ordering(C::AbstractMatrix{<:Integer}) -> Vector{Tuple{Symbol, Int}}, Vector{Int}

Given a Cartan matrix `C`, return a vector of `(family, rank)` pairs describing
its Dynkin type, together with a permutation vector `ord` giving a canonical
ordering of the rows/columns.

The algorithm decomposes the matrix into connected components via adjacency,
then classifies each component by its graph structure (path vs branching,
edge multiplicities).

Follows the same approach as OSCAR.jl's `cartan_type_with_ordering`.

# Examples
```julia
julia> cartan_type_with_ordering([2 -1; -1 2])
([(:A, 2)], [1, 2])

julia> cartan_type_with_ordering([2 -1; -2 2])
([(:B, 2)], [1, 2])
```
"""
function cartan_type_with_ordering(C::AbstractMatrix{<:Integer})
  rk = size(C, 1)
  @assert size(C, 1) == size(C, 2) "Cartan matrix must be square"

  type = Tuple{Symbol,Int}[]
  ord = sizehint!(Int[], rk)

  # Build adjacency list
  adj = [[j for j in 1:rk if i != j && C[i, j] != 0] for i in 1:rk]

  done = falses(rk)

  for v0 in 1:rk
    done[v0] && continue

    # ── Rank 1: isolated node ─────────────────────────────────────────
    if isempty(adj[v0])
      push!(type, (:A, 1))
      push!(ord, v0)
      done[v0] = true
      continue
    end

    # ── Rank 2: pair of nodes ─────────────────────────────────────────
    if length(adj[v0]) == 1 && length(adj[only(adj[v0])]) == 1
      v1 = only(adj[v0])
      prod = C[v0, v1] * C[v1, v0]
      if prod == 1
        push!(type, (:A, 2))
        push!(ord, v0, v1)
      elseif C[v0, v1] == -2
        # v0 is the short-root side → C_2 convention
        push!(type, (:C, 2))
        push!(ord, v0, v1)
      elseif C[v1, v0] == -2
        push!(type, (:B, 2))
        push!(ord, v0, v1)
      elseif C[v0, v1] == -3
        push!(type, (:G, 2))
        push!(ord, v0, v1)
      elseif C[v1, v0] == -3
        push!(type, (:G, 2))
        push!(ord, v1, v0)
      else
        error("Could not identify rank-2 Cartan matrix component")
      end
      done[v0] = true
      done[v1] = true
      continue
    end

    # ── Rank > 2: DFS to find the whole component ────────────────────
    comp = [v0]
    todo = [v0]
    done[v0] = true
    while !isempty(todo)
      v = pop!(todo)
      for w in adj[v]
        if !done[w]
          push!(comp, w)
          push!(todo, w)
          done[w] = true
        end
      end
    end
    sort!(comp)
    len_comp = length(comp)

    # Find degree-3 node (branching → D or E)
    deg3 = findfirst(v -> length(adj[v]) == 3, comp)

    if isnothing(deg3)
      # ── Path graph: A, B, C, or F ─────────────────────────────────
      # Find the start of the path (a leaf with simply-laced left neighbor)
      start = 0
      for v1 in filter(v -> length(adj[v]) == 1, comp)
        v2 = only(adj[v1])
        C[v1, v2] * C[v2, v1] == 1 || continue   # skip right end of B/C
        if len_comp == 4
          v3 = only(filter(!=(v1), adj[v2]))
          C[v2, v3] == -1 || continue               # skip right end of F
        end
        start = v1
        break
      end
      @assert start != 0 "Could not find start of path in component $comp"

      # Trace the path
      path = [start, only(adj[start])]
      for _ in 1:(len_comp - 2)
        push!(path, only(filter(!=(path[end - 1]), adj[path[end]])))
      end

      # Determine type from last edge
      if len_comp == 4 && C[path[3], path[2]] == -2
        push!(type, (:F, 4))
      elseif C[path[end - 1], path[end]] == -2
        push!(type, (:C, len_comp))
      elseif C[path[end], path[end - 1]] == -2
        push!(type, (:B, len_comp))
      else
        push!(type, (:A, len_comp))
      end
      append!(ord, path)
    else
      # ── Branching: D or E ──────────────────────────────────────────
      v_deg3 = comp[deg3]

      # Find the three paths from the branch node
      paths = [[v_deg3, v_n] for v_n in adj[v_deg3]]
      for path in paths
        while length(adj[path[end]]) == 2
          push!(path, only(filter(!=(path[end - 1]), adj[path[end]])))
        end
        popfirst!(path)  # remove the branch node itself
      end
      sort!(paths; by=length)

      @assert sum(length, paths) + 1 == len_comp

      if length(paths[2]) == 1
        # ── D type: two short arms of length 1 ──────────────────────
        push!(type, (:D, len_comp))
        if len_comp == 4
          push!(ord, only(paths[1]), v_deg3, only(paths[2]), only(paths[3]))
        else
          append!(ord, reverse!(paths[3]))
          push!(ord, v_deg3, only(paths[1]), only(paths[2]))
        end
      elseif length(paths[2]) == 2
        # ── E type: arms of length 1, 2, and 2/3/4 ─────────────────
        push!(type, (:E, len_comp))
        push!(ord, paths[2][2], only(paths[1]), paths[2][1], v_deg3)
        append!(ord, paths[3])
      else
        error("Could not identify branching Cartan matrix of rank $len_comp")
      end
    end
  end

  return type, ord
end

"""
    cartan_type(C::AbstractMatrix{<:Integer}) -> Vector{Tuple{Symbol, Int}}

Return the Cartan type of a Cartan matrix `C` as a vector of `(family, rank)` pairs.

# Examples
```julia
julia> cartan_type([2 -1; -1 2])
[(:A, 2)]

julia> cartan_type([2 -1 0 0; -1 2 0 0; 0 0 2 -1; 0 0 -2 2])
[(:A, 2), (:B, 2)]
```
"""
function cartan_type(C::AbstractMatrix{<:Integer})
  ct, _ = cartan_type_with_ordering(C)
  return ct
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Conversion: (Symbol, Int) pairs → Semisimple.jl DynkinType
# ═══════════════════════════════════════════════════════════════════════════════

"""
    _symbol_to_simple_type(fam::Symbol, rk::Int) -> Type{<:SimpleDynkinType}

Convert a `(:A, 3)` pair to `TypeA{3}`, etc.
"""
function _symbol_to_simple_type(fam::Symbol, rk::Int)
  fam === :A && return TypeA{rk}
  fam === :B && return TypeB{rk}
  fam === :C && return TypeC{rk}
  fam === :D && return TypeD{rk}
  fam === :E && return TypeE{rk}
  fam === :F && return TypeF4
  fam === :G && return TypeG2
  error("Unknown Dynkin family: $fam")
end

"""
    _cartan_type_to_dynkin_type(ct::Vector{Tuple{Symbol, Int}}) -> Type{<:DynkinType} or nothing

Convert a vector of `(family, rank)` pairs to a Semisimple.jl `DynkinType`.

Returns a single `SimpleDynkinType` if the vector has one element,
or a `ProductDynkinType` if it has multiple elements.
Returns `nothing` if the vector is empty.
"""
function _cartan_type_to_dynkin_type(ct::Vector{Tuple{Symbol,Int}})
  isempty(ct) && return nothing
  types = [_symbol_to_simple_type(fam, rk) for (fam, rk) in ct]
  length(types) == 1 && return types[1]
  return ProductDynkinType{Tuple{types...}}
end

# ═══════════════════════════════════════════════════════════════════════════════
#  String parsing: "A3", "B4", "A2xB3", etc. → DynkinType
# ═══════════════════════════════════════════════════════════════════════════════

"""
    parse_dynkin_type(s::AbstractString) -> Type{<:DynkinType}

Parse a string representation of a Dynkin type.

Components are separated by `x` or `×`. Each component is a letter
(`A`–`G`) followed by an integer rank.

# Examples
```julia
julia> parse_dynkin_type("A3")
TypeA{3}

julia> parse_dynkin_type("A2xB3")
ProductDynkinType{Tuple{TypeA{2}, TypeB{3}}}

julia> parse_dynkin_type("E6")
TypeE{6}
```
"""
function parse_dynkin_type(s::AbstractString)
  s = strip(String(s))
  isempty(s) && throw(ArgumentError("Empty Dynkin type string"))

  # Split on 'x' or '×'
  parts = split(s, r"[x×]")

  types = []
  for part in parts
    part = strip(part)
    isempty(part) && continue

    m = match(r"^([A-Ga-g])(\d+)$", part)
    m === nothing && throw(
      ArgumentError(
        "Cannot parse Dynkin type component: \"$part\". " *
        "Expected format like \"A3\", \"B4\", \"E6\".",
      ),
    )

    fam = Symbol(uppercase(m.captures[1]))
    rk = parse(Int, m.captures[2])

    # Validate
    _validate_cartan_type(fam, rk)

    push!(types, _symbol_to_simple_type(fam, rk))
  end

  isempty(types) && throw(ArgumentError("No valid Dynkin type components found in \"$s\""))

  if length(types) == 1
    return types[1]
  else
    return ProductDynkinType{Tuple{types...}}
  end
end

"""
    _validate_cartan_type(fam::Symbol, rk::Int)

Check that `(fam, rk)` is a valid Cartan type. Throws `ArgumentError` if not.
"""
function _validate_cartan_type(fam::Symbol, rk::Int)
  fam === :A && rk >= 1 && return nothing
  fam === :B && rk >= 2 && return nothing
  fam === :C && rk >= 2 && return nothing
  fam === :D && rk >= 4 && return nothing
  fam === :E && rk in (6, 7, 8) && return nothing
  fam === :F && rk == 4 && return nothing
  fam === :G && rk == 2 && return nothing
  throw(ArgumentError("Invalid Cartan type: ($fam, $rk)"))
end
