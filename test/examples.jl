@testset "Literate examples" begin
    EXAMPLEDIR = joinpath(@__DIR__, "..", "examples")
    for file in filter!(r"\.jl$", readdir(EXAMPLEDIR))
        @testset "$file" begin
            include(joinpath(EXAMPLEDIR, file))
        end
    end
end
