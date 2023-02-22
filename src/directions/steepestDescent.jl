struct SteepestDescent{F} <: ProgradioDirection{F}
    
    function SteepestDescent(; float_type::Type{F}=Float64) where {F<:AbstractFloat}
        return new{F}()
    end
end

mutable struct SteepestDescentState{F} <: ProgradioDirectionState{F}
    d::Vector{F}

    function SteepestDescentState(d::Vector{F};) where {F<:AbstractFloat}
        return new{F}(d)
    end
end

function direction_state(problem::ProgradioProblem{F}, ::SteepestDescent{F}) where {F<:AbstractFloat}
    return SteepestDescentState(zeros(F, length(problem.x_0)))
end

function direction!(state::IteratorState{F, DS, SS}, ::ProgradioProblem{F}, ::SteepestDescent{F}) where
    {F<:AbstractFloat, DS<:SteepestDescentState{F}, SS<:ProgradioSearchState}

    @. state.direction_state.d = -state.gx;
    return nothing
end

# Pretty printing
Base.show(io::IO, direction::SteepestDescent) = print(io, typeof(direction));