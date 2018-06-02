@testset "Literate examples" begin
    # We shouldn't run the examples that require Optim in Travis/CI,
    # because an update in LineSearches may be breaking with the
    # most recently tagged Optim version.
    if get(ENV, "CI", "") == "true"
        SKIPFILE = ["optim_linesearch.jl", "optim_initialstep.jl"]
    else
        SKIPFILE = ["",]
    end

    EXAMPLEDIR = joinpath(@__DIR__, "../docs/src/examples")

    myfilter(str) = r"\.jl$"(str) && !(str in SKIPFILE)
    for file in filter!(myfilter, readdir(EXAMPLEDIR))
        @testset "$file" begin
            mktempdir() do dir
                cd(dir) do
                    include(joinpath(EXAMPLEDIR, file))
                end
            end
        end
    end
end
