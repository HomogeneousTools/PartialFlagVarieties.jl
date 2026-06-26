#!/usr/bin/env julia
# ─────────────────────────────────────────────────────────────────────────────
#  Single source of truth for code formatting (JuliaFormatter, Blue style).
#
#  Run from anywhere in the repository:
#    julia devtools/format.jl           # format all tracked *.jl files in place
#    julia devtools/format.jl --check   # exit 1 if any tracked *.jl is unformatted
#
#  Both .github/workflows/format.yml (CI) and .hooks/pre-commit (local hook) call
#  this script, so local and CI formatting are byte-for-byte identical. To change
#  the formatter version, edit JULIAFORMATTER_VERSION below and NOTHING else —
#  pinning it stops a new upstream release from silently reformatting untouched
#  files and breaking CI.
# ─────────────────────────────────────────────────────────────────────────────

const JULIAFORMATTER_VERSION = "2.9.3"

using Pkg
# Install into a throwaway environment (explicit version => immune to any stale
# Manifest) so the project's own environment is left untouched.
Pkg.activate(; temp=true)
Pkg.add(
  Pkg.PackageSpec(; name="JuliaFormatter", version=JULIAFORMATTER_VERSION); io=devnull
)
using JuliaFormatter

# Operate from the repository root on git-tracked files only, so local scratch
# files are ignored and local/CI runs always see the same set.
cd(readchomp(`git rev-parse --show-toplevel`))
files = filter(f -> endswith(f, ".jl"), readlines(`git ls-files`))

check = "--check" in ARGS
# A comprehension (not `all(...)`) so every file is visited even in format mode.
results = [format(f; overwrite=(!check)) for f in files]

if check
  if all(results)
    println("✓ All ", length(files), " tracked Julia files are correctly formatted.")
  else
    bad = files[.!results]
    printstyled("JuliaFormatter check failed for ", length(bad), " file(s):\n"; color=:red)
    for f in bad
      println("  ", f)
    end
    println()
    println("Fix with:  julia devtools/format.jl  &&  git add -u")
    println("Skip with: git commit --no-verify")
    exit(1)
  end
end
