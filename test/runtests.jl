using Polyhedron 
using Test


@testset "funciton tests" begin

    @testset "is p-invariant" begin
        include("method_test.jl")
    end

    @testset "types" begin
        include("type_test.jl")
    end

end


