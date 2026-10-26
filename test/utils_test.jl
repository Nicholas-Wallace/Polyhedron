@testset "number of vectors ok" begin

    @test length(vet_eq_spc(1)) == 1
    @test length(vet_eq_spc(2)) == 2
    @test length(vet_eq_spc(3)) == 3

    @test length(vet_eq_spc(10)) == 10
    @test length(vet_eq_spc(11)) == 11

    @test_throws ErrorException vet_eq_spc(0)
    @test_throws ErrorException vet_eq_spc(-1)


end