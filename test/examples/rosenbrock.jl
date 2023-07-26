@testset "Box Rosenbrock" begin
    
    a = 1.0;
    b = 100.0;
    f(x) = (a - x[1])^2 + b * (x[2] - x[1]^2)^2;
    function g!(gx, x)
        gx[1] = 2 * (x[1] - a) + 4 * b * x[1] * (x[1]^2 - x[2]);
        gx[2] = 2 * b * (x[2] - x[1]^2);
        return nothing
    end

    box = Box([-2.0, -1.0], [0.0, 1.0]);
    x_guess = [-2.0, 1.0];

    nlp = NLProblem(x_guess, box, f, g!);

    sd = SteepestDescent();
    sd_armijo = TwoMetric(sd, Armijo());
    sd_armijo_solution = solve(nlp, sd_armijo);
    @test sd_armijo_solution.fx ≈ 1.0;
    sd_wolfe = TwoMetric(sd, Wolfe());
    sd_wolfe_solution = solve(nlp, sd_wolfe);
    @test sd_wolfe_solution.fx ≈ 1.0;

    fr = ConjugateGradient(10, 0.2, 10.0, FletcherReeves());
    fr_armijo = TwoMetric(fr, Armijo());
    fr_armijo_solution = solve(nlp, fr_armijo);
    @test fr_armijo_solution.fx ≈ 1.0;
    fr_wolfe = TwoMetric(fr, Wolfe());
    fr_wolfe_solution = solve(nlp, fr_wolfe);
    @test fr_wolfe_solution.fx ≈ 1.0;

    pr = ConjugateGradient(10, 0.2, 10.0, PolakRibiere());
    pr_armijo = TwoMetric(pr, Armijo());
    pr_armijo_solution = solve(nlp, pr_armijo);
    @test pr_armijo_solution.fx ≈ 1.0;
    pr_wolfe = TwoMetric(pr, Wolfe());
    pr_wolfe_solution = solve(nlp, pr_wolfe);
    @test pr_wolfe_solution.fx ≈ 1.0;

    hz = ConjugateGradient(10, 0.2, 10.0, HagerZhang());
    hz_armijo = TwoMetric(hz, Armijo());
    hz_armijo_solution = solve(nlp, hz_armijo);
    @test hz_armijo_solution.fx ≈ 1.0;
    hz_wolfe = TwoMetric(hz, Wolfe());
    hz_wolfe_solution = solve(nlp, hz_wolfe);
    @test hz_wolfe_solution.fx ≈ 1.0;

end