@testset "constraints added" begin

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
