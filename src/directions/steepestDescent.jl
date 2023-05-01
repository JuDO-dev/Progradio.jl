struct SteepestDescent <: Direction end

mutable struct SteepestDescentState{X<:AbstractVector} <: DirectionState
    d::X
end

build_state(x::AbstractVector, ::SteepestDescent) = SteepestDescentState(zero(x));

function direction!(state::IteratorState, ::SteepestDescent)
    
    @. state.algorithm_state.direction_state.d = -state.gx;
    
    return nothing
end