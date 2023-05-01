#= Empty problem, direction, search
struct EmptyProblem{F} <: P.ProgradioProblem{F} end
struct EmptyDirection{F} <: P.ProgradioDirection{F} end
struct EmptySearch{F} <: P.ProgradioSearch{F} end
empty_problem = EmptyProblem{Float}();
empty_direction = EmptyDirection{Float}();
empty_search = EmptySearch{Float}();

@testset "Iterator" begin
    
    i_max = 1000;
    x_tol = 1e-6;
    f_tol = 1e-7;
    g_tol = 1e-8;

    iterator = Iterator(empty_problem, empty_direction, empty_search;
        i_max=i_max, x_tol=x_tol, f_tol=f_tol, g_tol=g_tol);
    @test iterator.i_max == i_max;
    @test iterator.x_tol == x_tol;
    @test iterator.f_tol == f_tol;
    @test iterator.g_tol == g_tol;

end

# Empty states for problem, direction, search
struct EmptyDS{F} <: P.ProgradioDirectionState{F} end
struct EmptySS{F} <: P.ProgradioSearchState{F} end
empty_direction_state = EmptyDS{Float}();
empty_search_state = EmptySS{Float}();
n_x = 3;

empty_state = P.IteratorState(P.Iterating(), 0,
    zeros(n_x), Inf, zeros(n_x),
    zeros(n_x), 0.0, zeros(n_x),
    empty_direction_state, empty_search_state,
    falses(n_x), falses(n_x), falses(n_x), trues(n_x),
    0, 0.0
);

@testset "Exhaustion" begin
    
    i_max = 10;
    tol = 1e-6;
    iterator = Iterator(empty_problem, empty_direction, empty_search;
        i_max=i_max, x_tol=tol, f_tol=tol, g_tol=tol);
    empty_state.i = 10;

    next = iterate(iterator, empty_state);
    @test next isa Nothing;
    @test empty_state.status isa P.Exhausted;

end

@testset "Convergence" begin
    
    x_tol = 1e-6;
    empty_state.x_previous = ones(3);
    empty_state.x = ones(3);
    @test P.x_converged!(empty_state, x_tol, empty_problem);

    f_tol = 1e-6;
    empty_state.fx_previous = 1.0;
    empty_state.fx = 1.0;
    @test P.f_converged(empty_state, f_tol, empty_problem);

    g_tol = 1e-6;
    empty_state.gx = zeros(3);
    @test P.g_converged(empty_state, g_tol, empty_problem);

end=#