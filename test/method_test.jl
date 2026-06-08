@testset "Test 1" begin

    A = [1 1; 0 1]
    B = [2; 1]
    C = [1 0]
    X = [0.8 0; 0 1; -1 0; 0 -1]
    U = [1.2; -1.5]

    @test is_pinvariant(A, B, C, U, X) < 1

    @test is_pinvariant(A, B, C, U, X, SOF = true) > 1

    A = [1.2 0.2; -0.4 0.6]
    Ad = [-0.3 -0.2; 0.4 0.2]
    F = [-1 -1; 2 1]

    F = vcat(F, -F)
    @test is_pinvariant_delay(A, Ad, F, d = 1, symetric = false)["lambda"] < 1

end

@testset "Test symetric" begin

    A = [1.2 0.2; -0.4 0.6]
    Ad = [-0.3 -0.2; 0.4 0.2]
    F = [-1 -1; 2 1]

    @test is_pinvariant_delay(A, Ad, F, d = 1)["lambda"] < 1

end