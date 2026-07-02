# Contributing

## Use of LLMs

Parts of this package have been written with the assistance of large language
models. LLM-generated code is not inherently trustworthy: it can be subtly
wrong, miss edge cases, or introduce regressions. **Human review of every
change is essential** — please read and understand all code before merging,
regardless of how it was produced.

## Code formatting

The project uses [JuliaFormatter.jl](https://github.com/domluna/JuliaFormatter.jl)
with the **Blue** style (configured in `.JuliaFormatter.toml`).

> Formatting is driven by a single script,
> [`devtools/format.jl`](devtools/format.jl), which **pins** the JuliaFormatter
> version. CI and the pre-commit hook both call it, so local and CI formatting
> are byte-for-byte identical. **To change the formatter version, edit
> `JULIAFORMATTER_VERSION` in that script — the one and only place it lives.**
> (An unpinned install floats to the latest release and can silently reformat
> files you never touched, breaking CI — which is exactly what it once did.)

Format all tracked Julia files in place:

```bash
julia devtools/format.jl
```

Check formatting without modifying anything (exits non-zero if not):

```bash
julia devtools/format.jl --check
```

A **git pre-commit hook** is provided at `.hooks/pre-commit`. It runs
the check automatically before every commit that touches Julia files, so CI never
rejects your change due to formatting. Activate it once with:

```bash
git config core.hooksPath .hooks
```

If the hook fails, run the formatter, stage the result, and re-commit:

```bash
julia devtools/format.jl
git add -u
git commit ...
```

To bypass the hook on a work-in-progress commit use `git commit --no-verify`.
