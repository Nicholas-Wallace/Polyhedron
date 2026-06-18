using Test
using Polyhedron

@testset "is_pinvariant tests" begin
    
    A = [1 1; 0 1]
    B = [2; 1;;]
    C = [1 0]
    X = [0.8 0; 0 1; -1 0; 0 -1]
    U = [1.2; -1.5;;]
    
    # X é invariante por realimentação de estado
    @test is_pinvariant(A, B, C, U, X)["lambda"] < 1

    # X não é invariante por realimentação de saída
    @test_throws "Poliedro não é invariante" is_pinvariant(A, B, C, U, X, SOF=true)

    # Com Delay
    A = [1.2 0.2; -0.4 0.6]
    Ad = [-0.3 -0.2; 0.4 0.3]

    X = [-1 -1; 2 1]

    @test is_pinvariant(A, Ad, X, symetric=true, d=1)["lambda"] < 1
    @test is_pinvariant(A, Ad, X, symetric=true, d=10)["lambda"] < 1

    @test_throws "Poliedro não é invariante" @test is_pinvariant(A, Ad, X, symetric=true, d=11)["lambda"] < 1
end