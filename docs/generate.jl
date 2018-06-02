# generate examples
import Literate

EXAMPLEDIR = joinpath(@__DIR__, "src", "examples")
GENERATEDDIR = joinpath(@__DIR__, "src", "examples", "generated")
for example in filter!(r"\.jl$", readdir(EXAMPLEDIR))
    input = abspath(joinpath(EXAMPLEDIR, example))
    script = Literate.script(input, GENERATEDDIR)
    code = strip(read(script, String))
    mdpost(str) = replace(str, "@__CODE__" => code)
    Literate.markdown(input, GENERATEDDIR, postprocess = mdpost)
    Literate.notebook(input, GENERATEDDIR, execute = true)
end
