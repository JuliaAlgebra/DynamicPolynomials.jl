@testset "PolyVar order" begin
    @polyvar x
    z = x
    @polyvar x
    @test z != x
end
