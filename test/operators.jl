@testset "Projection" begin
    # Scalar
    @test P.Projection(0.0, -1.0, 1.0) == 0.0
    @test P.Projection(-2.0, -1.0, 1.0) == -1.0
    @test P.Projection(2.0, -1.0, 1.0) == 1.0
    @test P.Projection(1.0, -Inf, 0.0) == 0.0

    # Vector
    x_ℓ = [-1.0, 0.5, 0.0, -Inf];
    x_u = [-0.5, 1.0, 0.0, 0.0];
    @test P.Projection.(zeros(4), x_ℓ, x_u) == [-0.5, 0.5, 0.0, 0.0]
end

@testset "Binding" begin
    # Scalar
    @test P.Binding(-1.0, -1.0, 1.0, 0.1, 2.0) == true
    @test P.Binding(0.9, 0.0, 1.0, 0.1, -2.0) == true
    @test P.Binding(0.0, -1.0, 1.0, 0.1, 0.0) == false
    @test P.Binding(-1.0, -1.0, 1.0, 0.1, -2.0) == false
    
    # Vector
    x_ℓ = [0.0, 0.0, -1.0, -Inf];
    x_u = [0.0, 1.0, 1.0, 0.0];
    ϵ = [0.1, 0.1, 0.1, 0.1];
    ∇fx = [2.0, 2.0, 2.0, 2.0];
    @test P.Binding.(zeros(4), x_ℓ, x_u, ϵ, ∇fx) == [true, true, false, false]
end