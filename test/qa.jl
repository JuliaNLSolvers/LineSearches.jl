using LineSearches
using Test
import Aqua
import ExplicitImports
import JET

@testset "QA" begin
    @testset "Aqua" begin
        Aqua.test_all(LineSearches)
    end

    @testset "ExplicitImports" begin
        # No implicit imports (`using XY`)
        @test ExplicitImports.check_no_implicit_imports(LineSearches) === nothing

        # All explicit imports (`using XY: Z`) are loaded via their owners
        @test ExplicitImports.check_all_explicit_imports_via_owners(LineSearches) === nothing

        # No explicit imports (`using XY: Z`) that are not used
        @test ExplicitImports.check_no_stale_explicit_imports(LineSearches) === nothing

        # Nothing is accessed via modules other than its owner
        @test ExplicitImports.check_all_qualified_accesses_via_owners(LineSearches) === nothing

        # Almost no accesses of non-public names
        @test ExplicitImports.check_all_qualified_accesses_are_public(
            LineSearches;
            ignore = (
                :RefValue, # Base
                :max, # NaNMath
                :min, # NaNMath
            ),
        ) === nothing

        # No self-qualified accesses
        @test ExplicitImports.check_no_self_qualified_accesses(LineSearches) === nothing
    end

    @testset "JET" begin
        # Check that there are no undefined global references and undefined field accesses
        JET.test_package(LineSearches; target_defined_modules = true, mode = :typo)

        # Analyze methods based on their declared signature
        JET.test_package(LineSearches; target_defined_modules = true)
    end
end
