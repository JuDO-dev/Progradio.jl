@testset "Projection" begin
    # Scalar
    @test P.projection(0.0, -1.0, 1.0) == 0.0
    @test P.projection(-2.0, -1.0, 1.0) == -1.0
    @test P.projection(2.0, -1.0, 1.0) == 1.0
    @test P.projection(1.0, -Inf, 0.0) == 0.0

    # Vector
    x_ℓ = [-1.0, 0.5, 0.0, -Inf];
    x_u = [-0.5, 1.0, 0.0, 0.0];
    @test P.projection.(zeros(4), x_ℓ, x_u) == [-0.5, 0.5, 0.0, 0.0]
end

@testset "Binding" begin
    # Scalar
    @test P.binding(-1.0, -1.0, 1.0, 0.1, 2.0) == true
    @test P.binding(0.9, 0.0, 1.0, 0.1, -2.0) == true
    @test P.binding(0.0, -1.0, 1.0, 0.1, 0.0) == false
    @test P.binding(-1.0, -1.0, 1.0, 0.1, -2.0) == false
    
    # Vector
    x_ℓ = [0.0, 0.0, -1.0, -Inf];
    x_u = [0.0, 1.0, 1.0, 0.0];
    ϵ = [0.1, 0.1, 0.1, 0.1];
    gx = [2.0, 2.0, 2.0, 2.0];
    @test P.binding.(zeros(4), x_ℓ, x_u, ϵ, gx) == [true, true, false, false]
end