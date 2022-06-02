@testset "SBProblem type" begin
    # Example problem
    f(x) = sum(x.^2);
    function ∇f!(∇f, x)
        @. ∇f = 2 * x;
        return nothing
    end
    sbp = SBProblem(f, ∇f!, [1.5, -1.5], [1.0, -2.0], [2.0, -1.0]);

    # Tests
    @test sbp.x_0 == [1.5, -1.5]
    @test sbp.x_ℓ == [1.0, -2.0]
    @test sbp.x_u == [2.0, -1.0]
    @test sbp.f(sbp.x_0) == 4.5
    ∇f0 = zeros(2);
    sbp.∇f!(∇f0, sbp.x_0);
    @test ∇f0 == [3.0, -3.0]    
end

@testset "Quadratic example" begin
    qp = P.QP2([1.0, 0.5], [1.0, -0.5], [2.0, 0.5], 2.0);
    
    # Objective
    @test qp.f(qp.x_0) == 1.5
    @test qp.f(qp.x_ℓ) == 1.5
    @test qp.f(qp.x_u) == 4.5
    
    # Gradient
    ∇f0 = similar(qp.x_0);
    qp.∇f!(∇f0, qp.x_0);
    @test ∇f0 == [2.0, 2.0]
end

# @testset "Rosenbrock2" begin end