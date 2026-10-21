using Polyhedron
using Test

@testset "SLIT" begin
    n, m, p = 2, 1, 1
    A = [1.0 1.0; 0.0 1.0]
    B = reshape([2.0, 1.0], n, m)
    C = reshape([1.0, 0.0], p, n)
    G_ssf = [2; -2;;]
    G_sof = reshape([2.0], m, p)

    @testset "constructors" begin
        @testset "SLIT(A, B, C)  # defaults to OL" begin
            sys = SLIT(A, B, C)
            @test sys isa SLIT
            @test sys isa Polyhedron.System
            @test sys.A == A && sys.B == B && sys.C == C
            @test sys.feedback == OL
            @test size(sys.G) == (m, 0)
            @test sys.G == ones(m, 0)
        end

        @testset "SLIT(A, B, C, $fb)" for (fb, gshape) in
                ((SOF, (m, p)), (SSF, (m, n)), (OL, (m, 0)))
            sys = SLIT(A, B, C, fb)
            @test sys isa SLIT
            @test sys.feedback == fb
            @test sys.A == A && sys.B == B && sys.C == C
            @test size(sys.G) == gshape
            @test sys.G == ones(gshape...)           # default gain is all ones
        end

        @testset "initizalize with a G" begin
            sys = SLIT(A, B, C, SSF, G_ssf)
            @test sys isa SLIT
            @test sys.G == G_ssf

            sys = SLIT(A, B, C, SOF, G_sof)
            @test sys isa SLIT
            @test sys.G == G_sof
        end
    end

    @testset "dimension errors" begin
        @test_throws DimensionMismatch SLIT([1.0 0.0; 0.0 1.0; 0.0 0.0], B, C)  # A not square
        @test_throws DimensionMismatch SLIT(A, reshape([1.0, 0.0, 0.0], 3, 1), C) # B wrong rows
        @test_throws DimensionMismatch SLIT(A, B, reshape([1.0, 0.0, 0.0], 1, 3)) # C wrong cols
    end

    @testset "G is not validated (current behaviour)" begin
        # the constructor that receives a G is the auto generated one,
        # so the size of G is never checked
        sys = SLIT(A, B, C, SOF, G_ssf)     # SOF expects (m, p), G_ssf is (m, n)
        @test sys isa SLIT
        @test size(sys.G) == size(G_ssf)
        @test sys.G == G_ssf

        sys = SLIT(A, B, C, SSF, G_sof)     # SSF expects (m, n), G_sof is (m, p)
        @test sys isa SLIT
        @test size(sys.G) == size(G_sof)
        @test sys.G == G_sof

        # a Vector is not a Matrix{Float64}
        @test_throws MethodError SLIT(A, B, C, SOF, [2.0])

        # checkdims is skipped by this constructor as well
        @test SLIT([1.0 0.0; 0.0 1.0; 0.0 0.0], B, C, SSF, G_ssf) isa SLIT
    end

    @testset "other errors" begin
        @test_throws MethodError SLIT(A, B, C, "SOF")   # feedback must be a Feedback
        @test_throws MethodError SLIT(A, B)             # no 2-arg constructor
    end

    @testset "set_feedback!" begin
        sys = SLIT(A, B, C)                 # starts OL
        @test set_feedback!(sys, SOF) === sys
        @test sys.feedback == SOF
        @test size(sys.G) == (m, p)
        @test sys.G == ones(m, p)
        @test set_feedback!(sys, SSF) === sys
        @test sys.feedback == SSF
        @test size(sys.G) == (m, n)
        @test sys.G == ones(m, n)

        @testset "reset G" begin
            sys = SLIT(A, B, C, SSF, G_ssf)
            set_feedback!(sys, SOF)
            @test sys.G == ones(m, p)

            sys = SLIT(A, B, C, SOF, G_sof)
            set_feedback!(sys, SSF)
            @test sys.G == ones(m, n)
        end
    end
end

@testset "SystemDelayed" begin
    n, m, p = 2, 1, 1
    A = [1.0 1.0; 0.0 1.0]
    Ad = reshape([3.0, 4.0], n, m)      # input matrix (plays the role of B)
    C = reshape([1.0, 0.0], p, n)
    d = 2
    G_ssf = reshape([2.0, -2.0], m, n)
    G_sof = reshape([2.0], m, p)

    @testset "constructors" begin
        @testset "SystemDelayed(A, Ad, C, d)  # defaults to OL" begin
            sys = SystemDelayed(A, Ad, C, d)
            @test sys isa SystemDelayed
            @test sys isa Polyhedron.System
            @test sys.A == A && sys.Ad == Ad && sys.C == C
            @test sys.d == d
            @test sys.feedback == OL
            @test size(sys.G) == (m, 0)
            @test sys.G == ones(m, 0)
        end

        @testset "SystemDelayed(A, Ad, C, d, $fb)" for (fb, gshape) in
                ((SOF, (m, p)), (SSF, (m, n)), (OL, (m, 0)))
            sys = SystemDelayed(A, Ad, C, d, fb)
            @test sys isa SystemDelayed
            @test sys.feedback == fb
            @test sys.A == A && sys.Ad == Ad && sys.C == C
            @test sys.d == d
            @test size(sys.G) == gshape
            @test sys.G == ones(gshape...)            # default gain is all ones
        end

        @testset "number of inputs comes from Ad" begin
            Ad2 = reshape([1.0 2.0; 3.0 4.0], n, 2)   # two inputs
            @test size(SystemDelayed(A, Ad2, C, d, SOF).G) == (2, p)
            @test size(SystemDelayed(A, Ad2, C, d, SSF).G) == (2, n)
        end
    end

    @testset "dimension errors" begin
        @test_throws DimensionMismatch SystemDelayed([1.0 0.0; 0.0 1.0; 0.0 0.0], Ad, C, d)  # A not square
        @test_throws DimensionMismatch SystemDelayed(A, reshape([1.0, 0.0, 0.0], 3, 1), C, d) # Ad wrong rows
        @test_throws DimensionMismatch SystemDelayed(A, Ad, reshape([1.0, 0.0, 0.0], 1, 3), d) # C wrong cols
    end

    @testset "G is not validated (current behaviour)" begin
        sys = SystemDelayed(A, Ad, C, d, SOF, G_sof)
        @test sys isa SystemDelayed
        @test sys.G == G_sof

        sys = SystemDelayed(A, Ad, C, d, SOF, G_ssf)   # SOF expects (m, p), G_ssf is (m, n)
        @test sys isa SystemDelayed
        @test size(sys.G) == size(G_ssf)
        @test sys.G == G_ssf

        # a Vector is not a Matrix{Float64}
        @test_throws MethodError SystemDelayed(A, Ad, C, d, SOF, [2.0])

        # checkdims is skipped by this constructor as well
        @test SystemDelayed([1.0 0.0; 0.0 1.0; 0.0 0.0], Ad, C, d, SSF, G_ssf) isa SystemDelayed
    end

    @testset "other errors" begin
        @test_throws MethodError SystemDelayed(A, Ad, C, d, "SOF")  # feedback must be a Feedback
        @test_throws MethodError SystemDelayed(A, Ad, C)            # no 3-arg constructor
        @test_throws InexactError SystemDelayed(A, Ad, C, 1.5)      # d must be an Int
    end

    @testset "set_feedback!" begin
        sys = SystemDelayed(A, Ad, C, d)        # starts OL
        @test set_feedback!(sys, SOF) === sys
        @test sys.feedback == SOF
        @test size(sys.G) == (m, p)
        @test sys.G == ones(m, p)
        @test set_feedback!(sys, SSF) === sys
        @test sys.feedback == SSF
        @test size(sys.G) == (m, n)
        @test sys.G == ones(m, n)

        @testset "reset G" begin
            sys = SystemDelayed(A, Ad, C, d, SSF, G_ssf)
            set_feedback!(sys, SOF)
            @test sys.G == ones(m, p)

            sys = SystemDelayed(A, Ad, C, d, SOF, G_sof)
            set_feedback!(sys, SSF)
            @test sys.G == ones(m, n)
        end
    end
end

@testset "LTIDisturbance" begin
    n, m, p, q = 2, 1, 1, 1
    A = [1.0 1.0; 0.0 1.0]
    B = reshape([2.0, 1.0], n, m)
    E = reshape([0.5, -0.5], n, q)         # disturbance matrix
    C = reshape([1.0, 0.0], p, n)
    G_ssf = reshape([2.0, -2.0], m, n)
    G_sof = reshape([2.0], m, p)

    @testset "constructors" begin
        @testset "LTIDisturbance(A, B, E, C)  # defaults to OL" begin
            sys = LTIDisturbance(A, B, E, C)
            @test sys isa LTIDisturbance
            @test sys isa Polyhedron.System
            @test sys.A == A && sys.B == B && sys.E == E && sys.C == C
            @test sys.feedback == OL
            @test size(sys.G) == (m, 0)
            @test sys.G == ones(m, 0)
        end

        @testset "LTIDisturbance(A, B, E, C, $fb)" for (fb, gshape) in
                ((SOF, (m, p)), (SSF, (m, n)), (OL, (m, 0)))
            sys = LTIDisturbance(A, B, E, C, fb)
            @test sys isa LTIDisturbance
            @test sys.feedback == fb
            @test sys.A == A && sys.B == B && sys.E == E && sys.C == C
            @test size(sys.G) == gshape
            @test sys.G == ones(gshape...)            # default gain is all ones
        end

        @testset "E does not change the gain dimensions" begin
            E3 = reshape([0.5, -0.5, 1.0, 0.25, -0.25, 2.0], n, 3)   # three disturbances
            @test size(LTIDisturbance(A, B, E3, C, SOF).G) == (m, p)
            @test size(LTIDisturbance(A, B, E3, C, SSF).G) == (m, n)
        end
    end

    @testset "dimension errors" begin
        @test_throws DimensionMismatch LTIDisturbance([1.0 0.0; 0.0 1.0; 0.0 0.0], B, E, C)  # A not square
        @test_throws DimensionMismatch LTIDisturbance(A, reshape([1.0, 0.0, 0.0], 3, 1), E, C) # B wrong rows
        @test_throws DimensionMismatch LTIDisturbance(A, B, reshape([1.0, 0.0, 0.0], 3, 1), C) # E wrong rows
        @test_throws DimensionMismatch LTIDisturbance(A, B, E, reshape([1.0, 0.0, 0.0], 1, 3)) # C wrong cols
    end

    @testset "G is not validated (current behaviour)" begin
        sys = LTIDisturbance(A, B, E, C, SOF, G_sof)
        @test sys isa LTIDisturbance
        @test sys.G == G_sof

        sys = LTIDisturbance(A, B, E, C, SOF, G_ssf)   # SOF expects (m, p), G_ssf is (m, n)
        @test sys isa LTIDisturbance
        @test size(sys.G) == size(G_ssf)
        @test sys.G == G_ssf

        # a Vector is not a Matrix
        @test_throws MethodError LTIDisturbance(A, B, E, C, SOF, [2.0])

        # checkdims is skipped by this constructor as well
        @test LTIDisturbance([1.0 0.0; 0.0 1.0; 0.0 0.0], B, E, C, SSF, G_ssf) isa LTIDisturbance
    end

    @testset "other errors" begin
        @test_throws MethodError LTIDisturbance(A, B, E, C, "SOF")  # feedback must be a Feedback
        @test_throws MethodError LTIDisturbance(A, B, C)           # no 3-arg constructor
    end

    @testset "element type is not normalized" begin
        # the fields are ::Matrix, not ::Matrix{Float64}
        sys = LTIDisturbance([1 0; 0 1], [1; 0;;], [1; 0;;], [1 0])
        @test eltype(sys.A) == Int
        @test eltype(sys.G) == Float64
    end

    @testset "set_feedback!" begin
        sys = LTIDisturbance(A, B, E, C)      # starts OL
        @test set_feedback!(sys, SOF) === sys
        @test sys.feedback == SOF
        @test size(sys.G) == (m, p)
        @test sys.G == ones(m, p)
        @test set_feedback!(sys, SSF) === sys
        @test sys.feedback == SSF
        @test size(sys.G) == (m, n)
        @test sys.G == ones(m, n)

        @testset "reset G" begin
            sys = LTIDisturbance(A, B, E, C, SSF, G_ssf)
            set_feedback!(sys, SOF)
            @test sys.G == ones(m, p)

            sys = LTIDisturbance(A, B, E, C, SOF, G_sof)
            set_feedback!(sys, SSF)
            @test sys.G == ones(m, n)
        end
    end
end

@testset "LTIDisturbanceDelayedInput" begin
    n, m, p, q = 2, 1, 1, 1
    A = [1.0 1.0; 0.0 1.0]
    B = reshape([2.0, 1.0], n, m)
    E = reshape([0.5, -0.5], n, q)         # disturbance matrix
    C = reshape([1.0, 0.0], p, n)
    d = 2
    G_ssf = reshape([2.0, -2.0], m, n)
    G_sof = reshape([2.0], m, p)

    @testset "constructors" begin
        @testset "LTIDisturbanceDelayedInput(A, B, E, C, d)  # defaults to OL" begin
            sys = LTIDisturbanceDelayedInput(A, B, E, C, d)
            @test sys isa LTIDisturbanceDelayedInput
            @test sys isa Polyhedron.System
            @test sys.A == A && sys.B == B && sys.E == E && sys.C == C
            @test sys.d == d
            @test sys.feedback == OL
            @test size(sys.G) == (m, 0)
            @test sys.G == ones(m, 0)
        end

        @testset "LTIDisturbanceDelayedInput(A, B, E, C, d, $fb)" for (fb, gshape) in
                ((SOF, (m, p)), (SSF, (m, n)), (OL, (m, 0)))
            sys = LTIDisturbanceDelayedInput(A, B, E, C, d, fb)
            @test sys isa LTIDisturbanceDelayedInput
            @test sys.feedback == fb
            @test sys.A == A && sys.B == B && sys.E == E && sys.C == C
            @test sys.d == d
            @test size(sys.G) == gshape
            @test sys.G == ones(gshape...)            # default gain is all ones
        end

        @testset "E does not change the gain dimensions" begin
            E3 = reshape([0.5, -0.5, 1.0, 0.25, -0.25, 2.0], n, 3)   # three disturbances
            @test size(LTIDisturbanceDelayedInput(A, B, E3, C, d, SOF).G) == (m, p)
            @test size(LTIDisturbanceDelayedInput(A, B, E3, C, d, SSF).G) == (m, n)
        end
    end

    @testset "dimension errors" begin
        @test_throws DimensionMismatch LTIDisturbanceDelayedInput([1.0 0.0; 0.0 1.0; 0.0 0.0], B, E, C, d)  # A not square
        @test_throws DimensionMismatch LTIDisturbanceDelayedInput(A, reshape([1.0, 0.0, 0.0], 3, 1), E, C, d) # B wrong rows
        @test_throws DimensionMismatch LTIDisturbanceDelayedInput(A, B, reshape([1.0, 0.0, 0.0], 3, 1), C, d) # E wrong rows
        @test_throws DimensionMismatch LTIDisturbanceDelayedInput(A, B, E, reshape([1.0, 0.0, 0.0], 1, 3), d) # C wrong cols
    end

    @testset "G is not validated (current behaviour)" begin
        sys = LTIDisturbanceDelayedInput(A, B, E, C, d, SOF, G_sof)
        @test sys isa LTIDisturbanceDelayedInput
        @test sys.G == G_sof

        sys = LTIDisturbanceDelayedInput(A, B, E, C, d, SOF, G_ssf)   # SOF expects (m, p), G_ssf is (m, n)
        @test sys isa LTIDisturbanceDelayedInput
        @test size(sys.G) == size(G_ssf)
        @test sys.G == G_ssf

        # a Vector is not a Matrix
        @test_throws MethodError LTIDisturbanceDelayedInput(A, B, E, C, d, SOF, [2.0])

        # checkdims is skipped by this constructor as well
        @test LTIDisturbanceDelayedInput([1.0 0.0; 0.0 1.0; 0.0 0.0], B, E, C, d, SSF, G_ssf) isa LTIDisturbanceDelayedInput
    end

    @testset "other errors" begin
        @test_throws MethodError LTIDisturbanceDelayedInput(A, B, E, C, d, "SOF")  # feedback must be a Feedback
        @test_throws MethodError LTIDisturbanceDelayedInput(A, B, C, d)           # no 4-arg constructor
        @test_throws InexactError LTIDisturbanceDelayedInput(A, B, E, C, 1.5)    # d must be an Int
    end

    @testset "element type is not normalized" begin
        # the fields are ::Matrix, not ::Matrix{Float64}
        sys = LTIDisturbanceDelayedInput([1 0; 0 1], [1; 0;;], [1; 0;;], [1 0], 1)
        @test eltype(sys.A) == Int
        @test eltype(sys.G) == Float64
    end

    @testset "set_feedback!" begin
        sys = LTIDisturbanceDelayedInput(A, B, E, C, d)    # starts OL
        @test set_feedback!(sys, SOF) === sys
        @test sys.feedback == SOF
        @test size(sys.G) == (m, p)
        @test sys.G == ones(m, p)
        @test set_feedback!(sys, SSF) === sys
        @test sys.feedback == SSF
        @test size(sys.G) == (m, n)
        @test sys.G == ones(m, n)

        @testset "reset G" begin
            sys = LTIDisturbanceDelayedInput(A, B, E, C, d, SSF, G_ssf)
            set_feedback!(sys, SOF)
            @test sys.G == ones(m, p)

            sys = LTIDisturbanceDelayedInput(A, B, E, C, d, SOF, G_sof)
            set_feedback!(sys, SSF)
            @test sys.G == ones(m, n)
        end
    end
end

@testset "helpers" begin
    n, m, p, q = 2, 1, 1, 1
    A = [1.0 1.0; 0.0 1.0]
    B = reshape([2.0, 1.0], n, m)
    E = reshape([0.5, -0.5], n, q)
    C = reshape([1.0, 0.0], p, n)

    @testset "Feedback" begin
        @test length(instances(Feedback)) == 3
        @test Set(instances(Feedback)) == Set([SOF, SSF, OL])
        @test SOF != SSF && SSF != OL && SOF != OL
    end

    @testset "abstract System" begin
        @test isabstracttype(Polyhedron.System)
        @test SLIT <: Polyhedron.System
        @test SystemDelayed <: Polyhedron.System
        @test LTIDisturbance <: Polyhedron.System
        @test LTIDisturbanceDelayedInput <: Polyhedron.System
    end

    @testset "checkdims" begin
        @test checkdims(A, B, C) == (n, m, p)
        @test checkdims(A, B, E, C) == (n, m, p, q)
        @test_throws DimensionMismatch checkdims([1.0 0.0; 0.0 1.0; 0.0 0.0], B, C)      # A not square
        @test_throws DimensionMismatch checkdims(A, reshape([1.0, 0.0, 0.0], 3, 1), C)  # B wrong rows
        @test_throws DimensionMismatch checkdims(A, B, reshape([1.0, 0.0, 0.0], 1, 3))  # C wrong cols
        @test_throws DimensionMismatch checkdims(A, B, reshape([1.0, 0.0, 0.0], 3, 1), C)  # E wrong rows
    end

    @testset "gdims" begin
        @test Polyhedron.gdims(SOF, A, B, C) == (m, p)
        @test Polyhedron.gdims(SSF, A, B, C) == (m, n)
        @test Polyhedron.gdims(OL, A, B, C) == (m, 0)
    end
end
