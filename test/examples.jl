@testset "Literate examples" begin
    # We have to remove the Optim tests (and remove Optim from test/REQUIRE)
    # whenever LineSearches introduces a breaking change that the current
    # Optim release cannot handle.
    #   When the current Optim release works we should add the tests back.
    #SKIPFILE = ["optim_linesearch.jl", "optim_initialstep.jl"]
    SKIPFILE = []

    EXAMPLEDIR = joinpath(@__DIR__, "../docs/src/examples")
    myfilter(str) = occursin(r"\.jl$", str) && !(str in SKIPFILE)
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
