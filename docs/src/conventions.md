# Typesetting Conventions

This page documents the LaTeX conventions used throughout the
documentation and docstrings.

## Lie types and Dynkin diagrams

Lie type labels are set in upright roman type:
``\mathrm{A}_n``, ``\mathrm{B}_n``, ``\mathrm{C}_n``, ``\mathrm{D}_n``,
``\mathrm{E}_6``, ``\mathrm{E}_7``, ``\mathrm{E}_8``,
``\mathrm{F}_4``, ``\mathrm{G}_2``.

**Rationale**: the letters A–G denote names, not variables, so they
are upright.  Variable subscripts (rank ``n``) remain italic.

In source: `\\mathrm{A}_n`, `\\mathrm{E}_8`, etc.
In docs Markdown: `\mathrm{A}_n`, `\mathrm{E}_8`, etc.

## Hodge numbers

Hodge numbers are typeset as ``\mathrm{h}^{p,q}`` (roman h, superscript
``p,q``), consistent with the standard algebraic-geometry convention.

Betti numbers use roman ``\mathrm{b}_k``, Euler characteristic uses
``\chi`` (italic Greek letter).

In source docstrings: `\\mathrm{h}^{p,q}`.
In docs Markdown: `\mathrm{h}^{p,q}`.

## Geometric objects and operations

| Symbol | LaTeX | Meaning |
|:-------|:------|:--------|
| ``\mathrm{T}_{G/P}`` | `\mathrm{T}_{G/P}` | Tangent bundle of ``G/P`` |
| ``\Omega^p_{G/P}`` | `\Omega^p_{G/P}` | Bundle of ``p``-forms (``\Omega`` is upright by tradition) |
| ``\mathrm{K}_X`` | `\mathrm{K}_X` | Canonical class / canonical bundle |
| ``\mathrm{H}^k`` | `\mathrm{H}^k` | Cohomology functor (roman H) |
| ``\mathrm{H}^*`` | `\mathrm{H}^*` | Total cohomology |

Lowercase ``h^k = \dim \mathrm{H}^k`` is the Betti number; Hodge
number ``\mathrm{h}^{p,q} = \dim H^q(\Omega^p)`` has roman h.

## Group-theoretic operations

Use `\operatorname` for multi-letter operator names that are not
standard TeX commands:

| Symbol | LaTeX |
|:-------|:------|
| ``\operatorname{Z}(L)`` | `\operatorname{Z}(L)` |
| ``\operatorname{Pic}(X)`` | `\operatorname{Pic}(X)` |
| ``\operatorname{Sym}^k`` | `\operatorname{Sym}^k`` |
| ``\operatorname{End}`` | `\operatorname{End}` |

For the center of a reductive group we write ``\operatorname{Z}(L)^\circ``
(roman Z via `\operatorname`, superscript `\circ` for the connected component).

## Sheaves and bundles

| Symbol | LaTeX | Meaning |
|:-------|:------|:--------|
| ``\mathcal{O}`` | `\mathcal{O}` | Structure sheaf |
| ``\mathcal{O}(k)`` | `\mathcal{O}(k)` | Line bundle of degree ``k`` |
| ``\mathcal{S}`` | `\mathcal{S}` | Tautological subbundle |
| ``\mathcal{Q}`` | `\mathcal{Q}`` | Tautological quotient bundle |
| ``E^\vee`` | `E^\vee` | Dual bundle |

## Julia docstring escaping

In Julia `"""` docstrings, LaTeX must use **double backslash**:
`\\mathrm`, `\\operatorname`, etc., because the string is parsed by
Julia first. In Documenter `.md` files, single backslash is correct.

| In `.md` file | In Julia docstring |
|:--------------|:-------------------|
| `\mathrm{h}^{p,q}` | `\\mathrm{h}^{p,q}` |
| `\operatorname{Z}(L)` | `\\operatorname{Z}(L)` |
| `\mathrm{T}_{G/P}` | `\\mathrm{T}_{G/P}` |
