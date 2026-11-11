using Test
using Polyhedron

@testset "constraints added is_pinvariant SLIT" begin

    A = [1 1; 0 1]
    B = [2; 1;;]
    C = [1 0]
    X = [0.8 0; 0 1; -1 0; 0 -1]
    U = [1.2; -1.5;;]
    system = SLIT(A, B, C, SOF)

    model = is_pinvariant(system, X; U=U, return_model=true)

    @test constraint_by_name(model, "contractiave_pi0") != nothing
    @test constraint_by_name(model, "contractiave_pi1") != nothing
    @test constraint_by_name(model, "u_included0") != nothing
    @test constraint_by_name(model, "u_included1") != nothing

end

@testset "constraints added is_pinvariant SystemDelayed" begin

  A = [1.2 0.2; -0.4 0.6]
  Ad = [-0.3 -0.2; 0.4 0.2]
  C = [1 1]
  system = SystemDelayed(A, Ad, C, 1, SSF)

  X = [-1 -1; 2 1]
  X = vcat(X, -X)

  model = is_pinvariant(system, X; return_model=true)

  @test constraint_by_name(model, "contractiave_pi_H") != nothing
  @test constraint_by_name(model, "contractiave_pi_L") != nothing
  @test constraint_by_name(model, "contractiave_pi_M") != nothing
  @test constraint_by_name(model, "contractiave_pi_N") != nothing
  @test constraint_by_name(model, "contractiave_pi") != nothing

end

@testset "constraints added is_pinvariant LTIDisturbanceDelayedInput" begin

  A = [1.2 0.2; -0.4 0.6]
  B = [1.0; 0.5;;]
  E = [0.1 0.0; 0.0 0.1]
  C = [1 1]
  system = LTIDisturbanceDelayedInput(A, B, E, C, 1, SSF)

  X = [-1 -1; 2 1]
  X = vcat(X, -X)

  # find=false requires size(D, 1) == size(X, 1) for the current
  # `contractiave_pi` formulation
  D = [1.0 0.0; -1.0 0.0; 0.0 1.0; 0.0 -1.0]

  model = is_pinvariant(system, X, D; return_model=true)

  @test constraint_by_name(model, "contractiave_pi_H") != nothing
  @test constraint_by_name(model, "contractiave_pi_L") != nothing
  @test constraint_by_name(model, "contractiave_pi_M") != nothing
  @test constraint_by_name(model, "contractiave_pi_N") != nothing
  @test constraint_by_name(model, "contractiave_pi_P") != nothing
  @test constraint_by_name(model, "contractiave_pi_Q") != nothing
  @test constraint_by_name(model, "contractiave_pi")   != nothing

end
