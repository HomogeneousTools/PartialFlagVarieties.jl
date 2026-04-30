using Documenter
using Lie
using PartialFlagVarieties

makedocs(
  sitename="PartialFlagVarieties.jl",
  authors="Pieter Belmans",
  remotes=nothing,
  format=Documenter.HTML(
    prettyurls=get(ENV, "CI", nothing) == "true",
    canonical="https://pbelmans.github.io/PartialFlagVarieties.jl",
    assets=String[],
  ),
  modules=[PartialFlagVarieties],
  warnonly=[:missing_docs],
  pages=[
    "Home" => "index.md",
    "Getting Started" => "getting_started.md",
    "Conventions & Notation" => "conventions.md",
    "Common Workflows" => "workflows.md",
    "Mathematical Background" => "math.md",
    "API Reference" => [
      "MarkedDynkinType" => "api/marked_dynkin_type.md",
      "PartialFlagVariety" => "api/partial_flag_variety.md",
      "IrrepLevi" => "api/irrep_levi.md",
      "CompletelyReducibleBundle" => "api/bundle.md",
      "FilteredBundle" => "api/filtered_bundle.md",
      "UniversalBundles" => "api/universal_bundles.md",
      "Cohomology" => "api/cohomology.md",
      "Hodge & Hochschild" => "api/hodge.md",
      "ZeroLoci" => "api/zero_loci.md",
      "Koszul Algebra" => "api/koszul.md",
      "Constructions" => "api/constructions.md",
      "Labels" => "api/labels.md",
      "ExceptionalCollections" => "api/exceptional_collections.md",
      "Cache Configuration" => "api/cache_config.md",
    ],
  ],
  doctest=true,
)

deploydocs(
  repo="github.com/pbelmans/PartialFlagVarieties.jl.git",
  devbranch="main",
)
