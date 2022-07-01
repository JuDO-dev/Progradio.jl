@testset "Steepest Descent" begin
    struct SDIteratorState{F} <: P.ProgradioIteratorState{F}
        g::Vector{F}
        d::Vector{F}
    end

    state = SDIteratorState([1.0, -2.0, 0.0], zeros(3));

    P.memorize!(state, nothing);
    P.direction!(state, SteepestDescent());
    @test state.d == [-1.0, 2.0, 0.0]
end