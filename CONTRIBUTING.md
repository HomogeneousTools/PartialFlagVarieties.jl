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

> **The JuliaFormatter version is pinned** (currently `2.9.3`) in every command
> below, in `.github/workflows/format.yml`, and in `.hooks/pre-commit`. This
> keeps local, hook, and CI formatting byte-for-byte identical. An *unpinned*
> install floats to the latest release, which can silently reformat files you
> never touched and break CI. **When you bump the version, change it in all
> three places.**

Format all Julia files in-place:

```bash
julia -e 'using Pkg; Pkg.activate(; temp=true); Pkg.add(Pkg.PackageSpec(; name="JuliaFormatter", version="2.9.3")); using JuliaFormatter; format(".")'
```

Check whether the code is already correctly formatted (exits non-zero if not):

```bash
julia -e '
  using Pkg
  Pkg.activate(; temp=true)
  Pkg.add(Pkg.PackageSpec(; name="JuliaFormatter", version="2.9.3"))
  using JuliaFormatter
  exit(format(".", overwrite=false) ? 0 : 1)
'
```

A **git pre-commit hook** is provided at `.hooks/pre-commit`. It runs
the check automatically before every commit that touches Julia files, so CI never
rejects your change due to formatting. Activate it once with:

```bash
git config core.hooksPath .hooks
```

If the hook fails, run the formatter, stage the result, and re-commit:

```bash
julia -e 'using Pkg; Pkg.activate(; temp=true); Pkg.add(Pkg.PackageSpec(; name="JuliaFormatter", version="2.9.3")); using JuliaFormatter; format(".")'
git add -u
git commit ...
```

To bypass the hook on a work-in-progress commit use `git commit --no-verify`.
