@testset "BCIterator" begin
    # Example problem
    f(x) = sum(x .^ 2);
    function g!(g, x)
        @. g = 2 * x;
        return nothing
    end
    bcp = BCProblem(f, g!, [1.0, -2.0], [2.0, -1.0], [1.5, -1.5], 0.01);

    struct EmptyDirection{F} <: P.ProgradioDirection{F} end
    struct EmptyOptimizer{F, D} <: P.ProgradioOptimizer{F, D} end
    optimizer = EmptyOptimizer{Float64, EmptyDirection{Float64}}();

    # iterator
    bci = iterator(bcp, optimizer);
    @test bci.i_max == 20
    @test bci.f_tol == 1e-6
    @test bci.g_tol == 1e-6
    @test bci.x_tol == 1e-6
end