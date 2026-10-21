using Test
using Polyhedron

@testset "is_pinvariant tests" begin
    
    @testset "LTI Systems" begin
    # sistema do book abstract
    A = [1 1; 0 1]
    B = [2; 1;;]
    C = [1 0]
    X = [0.8 0; 0 1; -1 0; 0 -1]
    U = [1.2; -1.5;;]

    @testset "SSF" begin
	system = SLIT(A, B, C, SSF)
    
        # X é invariante por realimentação de estado
        ip, matrices = is_pinvariant(system, U, X, verbose=true)

	@test ip.λ < 1
	@test ip.system == system
	
	# testar as condições
	lx = size(X,1); lu = size(U,1)
	xi = ones(lx); phi = ones(lu)
	@test matrices["H"] * matrices["X"] ≈ matrices["X"] * (system.A + system.B * matrices["F"]) atol=0.1
	@test matrices["M"] * matrices["X"] ≈ matrices["U"] * matrices["F"] atol=0.1
	@test all(matrices["H"] * xi .<= xi * matrices["lambda"])
	@test all(matrices["M"] * xi .<= phi)

	# TODO add a test for a not invariant by SSF polyhedron
    end

    @testset "SOF" begin

        # X não é invariante por realimentação de saída
	system = SLIT(A, B, C, SOF)
    	@test_throws "Poliedro não é invariante" is_pinvariant(system, U, X)
    end
    end

    @testset "SystemDelayed" begin
    # Com Delay
    A = [1.2 0.2; -0.4 0.6]
    Ad = [-0.3 -0.2; 0.4 0.3]

    X = [-1 -1; 2 1]

    # TODO testes com delay
    #@test is_pinvariant(SystemDelayed(A, Ad, 1), X, symetric=true)["lambda"] < 1
    #@test is_pinvariant(SystemDelayed(A, Ad, 10), X, symetric=true)["lambda"] < 1

    #@test_throws "Poliedro não é invariante" @test is_pinvariant(SystemDelayed(A, Ad, 11), X, symetric=true)["lambda"] < 1
    end
end
