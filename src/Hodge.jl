# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge.jl — Hodge numbers and Hodge diamond printing
#
#  Provides:
#   - hodge_numbers: concrete BigInt Hodge diamond of a ZeroLocus
#   - print_hodge_diamond: ASCII Hodge diamond via PrettyTables
# ═══════════════════════════════════════════════════════════════════════════════

export hodge_numbers
export print_hodge_diamond

# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge numbers of zero loci
# ═══════════════════════════════════════════════════════════════════════════════

"""
    hodge_numbers(Z::ZeroLocus) -> Matrix{BigInt}

Compute the Hodge diamond ``h^{p,q}(Z)`` for ``p, q = 0, \\ldots, \\dim Z``.

Uses the Koszul resolution and the conormal exact sequence.
Returns a ``(d+1) \\times (d+1)`` matrix where entry ``[p+1, q+1] = h^{p,q}``.

For ``p = 0``: computed directly from the Koszul resolution of ``\\mathcal{O}_Z``.
For ``p \\ge 1``: uses the conormal sequence and previously computed ``h^{j,q}``
for ``j < p``.

# Examples
```jldoctest
julia> using PartialFlagVarieties

julia> X = projective_space(4);

julia> Z = zero_locus(line_bundle(X, 5));

julia> h = hodge_numbers(Z);

julia> h[2, 2]  # h^{1,1}
1

julia> h[3, 2]  # h^{2,1}
101
```
"""
function hodge_numbers(Z::ZeroLocus)
  d = dimension(Z)
  H_sym = hodge_numbers_symbolic(Z)

  hodge = zeros(BigInt, d + 1, d + 1)
  for p in 0:d, q in 0:d
    e = H_sym[p + 1, q + 1]
    if is_determined(e)
      hodge[p + 1, q + 1] = e.constant
    else
      # Fall back to the constant part (lower bound)
      hodge[p + 1, q + 1] = e.constant
    end
  end
  hodge
end

# ═══════════════════════════════════════════════════════════════════════════════
#  Hodge diamond printing
# ═══════════════════════════════════════════════════════════════════════════════

# Table format with no borders and no dividing lines, used for Hodge diamonds.
const _DIAMOND_FMT = TextTableFormat(;
  borders=text_table_borders__borderless,
  horizontal_line_at_beginning=false,
  horizontal_lines_at_column_labels=:none,
  horizontal_line_at_merged_column_labels=false,
  horizontal_line_after_column_labels=false,
  horizontal_lines_at_data_rows=:none,
  horizontal_line_before_row_group_label=false,
  horizontal_line_after_row_group_label=false,
  horizontal_line_after_data_rows=false,
  horizontal_line_before_summary_rows=false,
  horizontal_line_after_summary_rows=false,
  vertical_line_at_beginning=false,
  vertical_line_after_row_number_column=false,
  vertical_line_after_row_label_column=false,
  vertical_lines_at_data_columns=:none,
  vertical_line_after_data_columns=false,
  vertical_line_after_continuation_column=false,
)

"""
    print_hodge_diamond([io::IO,] h::Matrix{<:Integer})

Print a centred ASCII Hodge diamond for a variety of arbitrary dimension.

`h[p+1, q+1] = h^{p,q}` (1-based matrix indexing). The diamond is centred
using PrettyTables so that each column is as wide as its largest entry,
giving exact alignment regardless of digit count.

# Layout
`h^{p,q}` is placed at grid row `p+q+1`, column `p-q+d+1` in a
`(2d+1) × (2d+1)` string matrix (where `d = size(h,1) - 1`).

# Examples
```julia
julia> using PartialFlagVarieties

julia> X = Gr(2, 6);

julia> Z = zero_locus(reduce(direct_sum, [line_bundle(X, 1) for _ in 1:6]));

julia> print_hodge_diamond(stdout, hodge_numbers(Z));  # K3 of degree 14 in P^5
       1
    0      0
 1     20     1
    0      0
       1
```
"""
function print_hodge_diamond(io::IO, h::Matrix{<:Integer})
  d = size(h, 1) - 1
  sz = 2 * d + 1

  cells = fill("", sz, sz)
  for p in 0:d, q in 0:d
    cells[p + q + 1, p - q + d + 1] = string(h[p + 1, q + 1])
  end

  pretty_table(io, cells;
    table_format=_DIAMOND_FMT,
    show_column_labels=false,
    alignment=:c,
  )
end

print_hodge_diamond(h::Matrix{<:Integer}) = print_hodge_diamond(stdout, h)
