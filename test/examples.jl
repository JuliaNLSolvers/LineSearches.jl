@testset "Literate examples" begin
    # We shouldn't run the examples that require Optim in Travis/CI,
    # because an update in LineSearches may be breaking with the
    # most recently tagged Optim version.
    noCI = ["optim_linesearch.jl", "optim_initialstep.jl"] # used in generate.jl

    EXAMPLEDIR = joinpath(@__DIR__, "../docs/src/examples")
    for file in filter!(r"\.jl$", readdir(EXAMPLEDIR))
        if get(ENV, "CI", "") == "true"
            file in noCI && continue
        end
        @testset "$file" begin
            mktempdir() do dir
                cd(dir) do
                    include(joinpath(EXAMPLEDIR, file))
                end
            end
        end
    end
end
