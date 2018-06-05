@testset "Literate examples" begin
    # TODO: Remove items from `SKIPFILE` as soon as they run on the latest
    # stable `Optim` (or other dependency)
    SKIPFILE = ["optim_linesearch.jl", "optim_initialstep.jl"]

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
