@testset "BCProblem" begin

    #Example QP
    x_0 = ones(3);
    ℓ = zeros(3);
    u = ones(3);
    f(x) = P.dot(x, x);
    function g!(gx, x)
        @. gx = 2 * x;
        return nothing
    end
    bcp = BCProblem(x_0, ℓ, u, f, g!);
    @test bcp.x_0 == x_0;
    @test bcp.ℓ == ℓ;
    @test bcp.u == u;
    @test bcp.f(bcp.x_0) == 3.0;
    gx_0 = zeros(3);
    bcp.g!(gx_0, bcp.x_0);
    @test gx_0 == 2.0 * ones(3);

    # Box dimensions match
    ℓ = zeros(4);
    @test_throws DimensionMismatch BCProblem(x_0, ℓ, u, f, g!)
    ℓ = zeros(3);
    u = ones(4);
    @test_throws DimensionMismatch BCProblem(x_0, ℓ, u, f, g!)

    # x_0 inside box
    ℓ = zeros(3);
    u = zeros(3);
    x_0 = ones(3);
    @test_throws DomainError BCProblem(x_0, ℓ, u, f, g!)

    # Non-negative B_tol
    @test_throws DomainError BCProblem(x_0, ℓ, u, -1.0, f, g!)

    # Binding tolerance 10% of smallest positive bound range
    x_0 = zeros(3);
    ℓ = zeros(3);
    u = [10.0, 1.0, 0.0];
    bcp = BCProblem(x_0, ℓ, u, f, g!);
    @test bcp.B_tol == 0.1;
    
end