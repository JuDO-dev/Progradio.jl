@testset "Steepest Descent" begin
    
    sd = SteepestDescent(float_type=Float, integer_type=Int);
    d = ones(Float, 3);
    sd_state = P.SteepestDescentState(d; integer_type=Int);
    gx = ones(Float, 3);
    
    test_state = P.IteratorState(P.Iterating(), 0,
        zeros(n_x), Inf, d,
        zeros(n_x), 0.0, zeros(n_x),
        sd_state, empty_search_state,
        falses(n_x), falses(n_x), falses(n_x), trues(n_x),
        0, 0.0
    );

    P.direction!(test_state, sd);
    @test test_state.direction_state.d == [-1.0, -1.0, -1.0];

end