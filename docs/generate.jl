# generate examples
import Literate

# We shouldn't run the examples that require Optim in Travis/CI,
# because an update in LineSearches may be breaking with the
# most recently tagged Optim version.
if get(ENV, "CI", "") == "true"
    ONLYSTATIC = ["optim_linesearch.jl", "optim_initialstep.jl"]
else
    ONLYSTATIC = ["",]
end

EXAMPLEDIR = joinpath(@__DIR__, "src", "examples")
GENERATEDDIR = joinpath(@__DIR__, "src", "examples", "generated")
for example in filter!(r"\.jl$", readdir(EXAMPLEDIR))
    input = abspath(joinpath(EXAMPLEDIR, example))
    script = Literate.script(input, GENERATEDDIR)
    code = strip(read(script, String))
    mdpost(str) = replace(str, "@__CODE__" => code)
    Literate.markdown(input, GENERATEDDIR, postprocess = mdpost,
                      documenter = !(example in ONLYSTATIC))
    Literate.notebook(input, GENERATEDDIR, execute = !(example in ONLYSTATIC))
end
