#=@testset "Steepest Descent" begin
    
    sd = SteepestDescent(float_type=Float);
    d = ones(Float, 3);
    sd_state = P.SteepestDescentState(d);
    gx = ones(Float, 3);
    
    test_state = P.IteratorState(P.Iterating(), 0,
        zeros(n_x), Inf, d,
        zeros(n_x), 0.0, zeros(n_x),
        sd_state, empty_search_state,
        falses(n_x), falses(n_x), falses(n_x), trues(n_x),
        0, 0.0
    );

    struct TestProblem{F} <: P.ProgradioProblem{F} end
    test_problem = TestProblem{Float}();

    P.direction!(test_state, test_problem, sd);
    @test test_state.direction_state.d == [-1.0, -1.0, -1.0];

end=#