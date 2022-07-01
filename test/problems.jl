@testset "BCProblem type" begin
    # Example problem
    f(x) = sum(x .^ 2);
    function g!(g, x)
        @. g = 2 * x;
        return nothing
    end
    bcp = BCProblem(f, g!, [1.0, -2.0], [2.0, -1.0], [1.5, -1.5], 0.01);

    # Tests
    @test bcp.f(bcp.x_0) == 4.5
    g0 = zeros(2);
    bcp.g!(g0, bcp.x_0);
    @test g0 == [3.0, -3.0]  
    @test bcp.x_ℓ == [1.0, -2.0]
    @test bcp.x_u == [2.0, -1.0]
    @test bcp.x_0 == [1.5, -1.5]
    @test bcp.ϵ == 0.01
end