@testset "Steepest Descent" begin
    sd = SteepestDescent();
    sdstate = P.direction_state(3, sd);
    
    struct TestOptimizerState{F, DS} <: P.ProgradioOptimizerState{F, DS}
        gx::Vector{F}
        d::Vector{F}
    end
    state = TestOptimizerState{Float64, P.SteepestDescentState{Float64}}([1.0, -2.0, 0.0], zeros(3));

    P.direction!(state, sd);
    @test state.d == [-1.0, 2.0, 0.0]
end