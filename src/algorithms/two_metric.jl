struct TwoMetric{D<:Direction, LS<:LineSearch} <: Algorithm
    direction::D
    line_search::LS
end

mutable struct TwoMetricState{DS<:DirectionState, LSS<:LineSearchState} <: AlgorithmState
    direction_state::DS
    line_search_state::LSS
end

build_state(x::AbstractArray, tm::TwoMetric) = TwoMetricState(build_state(x, tm.direction), build_state(x, tm.line_search));

function iterate_algorithm!(state::IteratorState, problem, algorithm)

    direction!(state, algorithm.direction);

    line_search!(state, problem, algorithm.line_search);

    problem.g!(state.gx, state.x);

    binding!(state, problem.constraints);

    return nothing
end