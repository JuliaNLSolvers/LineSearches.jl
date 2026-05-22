if Base.HOME_PROJECT[] !== nothing
    # JuliaLang/julia/pull/28625
    Base.HOME_PROJECT[] = abspath(Base.HOME_PROJECT[])
end

using Documenter, LineSearches

# Generate examples
include("generate.jl")

GENERATEDEXAMPLES = [joinpath("examples", "generated", f) for f in (
    "customoptimizer.md", "optim_linesearch.md", "optim_initialstep.md")]

# Build documentation.
makedocs(
    format = Documenter.HTML(),
    sitename = "LineSearches.jl",
    doctest = false,
    pages = Any[
        "Home" => "index.md",
        "Examples" => GENERATEDEXAMPLES,
        "API Reference" => [
            "reference/linesearch.md",
            "reference/initialstep.md",
            ]
        ]
    )

deploydocs(
    repo = "github.com/JuliaNLSolvers/LineSearches.jl.git",
)
