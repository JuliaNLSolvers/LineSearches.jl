@testset "Literate examples" begin
    EXAMPLEDIR = joinpath(@__DIR__, "../docs/src/examples")
    for file in filter!(r"\.jl$", readdir(EXAMPLEDIR))
        @testset "$file" begin
            mktempdir() do dir
                cd(dir) do
                    include(joinpath(EXAMPLEDIR, file))
                end
            end
        end
    end
end
