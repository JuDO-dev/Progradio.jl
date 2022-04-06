@testset "Projection" begin
    # Scalar
    @test Projection(0.0, -1.0, 1.0) == 0.0
    @test Projection(-2.0, -1.0, 1.0) == -1.0
    @test Projection(2.0, -1.0, 1.0) == 1.0
    @test Projection(1.0, -Inf, 0.0) == 0.0

    # Vector
    ℓ = [-1.0, 0.5, 0.0, -Inf];
    u = [-0.5, 1.0, 0.0, 0.0];
    @test Projection.(zeros(4), ℓ, u) == [-0.5, 0.5, 0.0, 0.0]
end

@testset "Binding" begin
    # Scalar
    @test Binding(-1.0, -1.0, 1.0, 0.1, 2.0) == true
    @test Binding(0.9, 0.0, 1.0, 0.1, -2.0) == true
    @test Binding(0.0, -1.0, 1.0, 0.1, 0.0) == false 
    @test Binding(-1.0, -1.0, 1.0, 0.1, -2.0) == false
    
    # Vector
    ℓ = [0.0, 0.0, -1.0, -Inf];
    u = [0.0, 1.0, 1.0, 0.0];
    ϵ = [0.1, 0.1, 0.1, 0.1];
    ∇fx = [2.0, 2.0, 2.0, 2.0];
    @test Binding.(zeros(4), ℓ, u, ϵ, ∇fx) == [true, true, false, false]
end