# generate examples
import Literate

# We shouldn't run the examples that require Optim in Travis/CI,
# because an update in LineSearches may be breaking with the
# most recently tagged Optim version.
noCI = ["optim_linesearch.jl", "optim_initialstep.jl"] # used in generate.jl

EXAMPLEDIR = joinpath(@__DIR__, "src", "examples")
GENERATEDDIR = joinpath(@__DIR__, "src", "examples", "generated")
for example in filter!(r"\.jl$", readdir(EXAMPLEDIR))
    if get(ENV, "CI", "") == "true"
        example in noCI && continue
    end
    input = abspath(joinpath(EXAMPLEDIR, example))
    script = Literate.script(input, GENERATEDDIR)
    code = strip(read(script, String))
    mdpost(str) = replace(str, "@__CODE__" => code)
    Literate.markdown(input, GENERATEDDIR, postprocess = mdpost)
    Literate.notebook(input, GENERATEDDIR, execute = true)
end
