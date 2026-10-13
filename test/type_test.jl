using Polyhedron
using Test

@testset "SLIT" begin
    n, m, p = 2, 1, 1
    A = [1.0 1.0; 0.0 1.0]
    B = reshape([2.0, 1.0], n, m)
    C = reshape([1.0, 0.0], p, n)

    @testset "constructors" begin
        @testset "SLIT(A, B, C)  # defaults to OL" begin
            sys = SLIT(A, B, C)
            @test sys isa SLIT
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
    end

    @testset "dimension errors" begin
        @test_throws DimensionMismatch SLIT([1.0 0.0; 0.0 1.0; 0.0 0.0], B, C)  # A not square
        @test_throws DimensionMismatch SLIT(A, reshape([1.0, 0.0, 0.0], 3, 1), C) # B wrong rows
        @test_throws DimensionMismatch SLIT(A, B, reshape([1.0, 0.0, 0.0], 1, 3)) # C wrong cols
    end

    @testset "other errors" begin
        @test_throws MethodError SLIT(A, B, C, "SOF")   # feedback must be a Feedback
        @test_throws MethodError SLIT(A, B)             # no 2-arg constructor
    end

    @testset "set_feedback!" begin
        sys = SLIT(A, B, C)                 # starts OL
        set_feedback!(sys, SOF)
        @test sys.feedback == SOF
        @test size(sys.G) == (m, p)
        set_feedback!(sys, SSF)
        @test sys.feedback == SSF
        @test size(sys.G) == (m, n)
    end
end
