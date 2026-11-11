@testset "constraints added SLIT" begin

    A = [1 1; 0 1]
    B = [2; 1;;]
    C = [1 0]
    X = [0.8 0; 0 1; -1 0; 0 -1]
    U = [1.2; -1.5;;]
    system = SLIT(A, B, C, SOF)

    model = pinvariant_poly(system, X; U=U, return_model=true)

    @test constraint_by_name(model, "contractiave_pi0") != nothing
    @test constraint_by_name(model, "contractiave_pi1") != nothing
    @test constraint_by_name(model, "x_included0") != nothing
    @test constraint_by_name(model, "x_included1") != nothing
    @test constraint_by_name(model, "u_included0") != nothing
    @test constraint_by_name(model, "u_included1") != nothing
    @test constraint_by_name(model, "pseudo_inv") != nothing

    for i in 1:8
        @test constraint_by_name(model, "enl[$i]") != nothing
    end

    model_enl = pinvariant_poly(system, X; U=U, enl=state_constraint_enl!)

    @test constraint_by_name(model_enl, "state_enl0") != nothing
    @test constraint_by_name(model_enl, "state_enl1") != nothing

end

@testset "constraints added SystemDelayed" begin

  A = [1.2 0.2; -0.4 0.6]
	Ad = [-0.3 -0.2; 0.4 0.2]
	C = [1 1]
  system = SystemDelayed(A, Ad, C, 1, SSF)

  X = [-1 -1; 2 1]
  X = vcat(X, -X)

  model = pinvariant_poly(system, X; return_model=true)

  @test constraint_by_name(model, "contractiave_pi_H") != nothing
  @test constraint_by_name(model, "contractiave_pi_L") != nothing
  @test constraint_by_name(model, "contractiave_pi_M") != nothing
  @test constraint_by_name(model, "contractiave_pi_N") != nothing
  @test constraint_by_name(model, "contractiave_pi") != nothing

  @test constraint_by_name(model, "x_included0") != nothing
  @test constraint_by_name(model, "x_included1") != nothing

  for i in 1:8
    @test constraint_by_name(model, "enl[$i]") != nothing
  end

  model_enl = pinvariant_poly(system, X; enl=state_constraint_enl!, return_model=true)

  @test constraint_by_name(model_enl, "state_enl0") != nothing
  @test constraint_by_name(model_enl, "state_enl1") != nothing

end

@testset "constraints added LTIDisturbanceDelayedInput" begin

  A = [1.2 0.2; -0.4 0.6]
  B = [1.0; 0.5;;]
  E = [0.1 0.0; 0.0 0.1]
  C = [1 1]
  system = LTIDisturbanceDelayedInput(A, B, E, C, 1, SSF)

  X = [-1 -1; 2 1]
  X = vcat(X, -X)

  D = [1.0 0.0; -1.0 0.0; 0.0 1.0; 0.0 -1.0; 1.0 1.0; -1.0 -1.0]

  model = pinvariant_poly(system, X, D; return_model=true)

  @test constraint_by_name(model, "contractiave_pi_H") != nothing
  @test constraint_by_name(model, "contractiave_pi_L") != nothing
  @test constraint_by_name(model, "contractiave_pi_M") != nothing
  @test constraint_by_name(model, "contractiave_pi_N") != nothing
  @test constraint_by_name(model, "contractiave_pi_P") != nothing
  @test constraint_by_name(model, "contractiave_pi_Q") != nothing
  @test constraint_by_name(model, "contractiave_pi")   != nothing

  @test constraint_by_name(model, "x_included0") != nothing
  @test constraint_by_name(model, "x_included1") != nothing

  for i in 1:8
    @test constraint_by_name(model, "enl[$i]") != nothing
  end

  model_enl = pinvariant_poly(system, X, D; enl=state_constraint_enl!, return_model=true)

  @test constraint_by_name(model_enl, "state_enl0") != nothing
  @test constraint_by_name(model_enl, "state_enl1") != nothing

end
