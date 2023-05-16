struct SteepestDescent <: Direction end

mutable struct SteepestDescentState{T} <: DirectionState{T}
    d::Vector{T}
end

build_state(x::X, ::SteepestDescent) where {T<:Real, N, X<:AbstractArray{T, N}} =
    SteepestDescentState(Vector{T}(undef, length(x))
);

function direction!(state::IteratorState, ::SteepestDescent)
    
    @. state.algorithm_state.direction_state.d = -state.gx;
    
    return nothing
end