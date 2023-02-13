struct SteepestDescent{F, I} <: ProgradioDirection{F, I}
    
    function SteepestDescent(; float_type::Type{F}=Float64, integer_type::Type{I}=Int64) where {F<:AbstractFloat, I<:Integer}
        return new{F, I}()
    end
end

mutable struct SteepestDescentState{F, I} <: ProgradioDirectionState{F, I}
    d::Vector{F}

    function SteepestDescentState(d::Vector{F}; integer_type::Type{I}=Int64) where {F<:AbstractFloat, I<:Integer}
        return new{F, I}(d)
    end
end

function direction_state(problem::ProgradioProblem{F, I}, ::SteepestDescent{F, I}) where {F<:AbstractFloat, I<:Integer}
    return SteepestDescentState(zeros(F, length(problem.x_0)), integer_type=I)
end

function direction!(state::IteratorState{F, I, DS, SS}, ::SteepestDescent{F, I}) where
    {F<:AbstractFloat, I<:Integer, DS<:SteepestDescentState{F, I}, SS<:ProgradioSearchState}

    @. state.direction_state.d = -state.gx;
    return nothing
end

# Pretty printing
Base.show(io::IO, direction::SteepestDescent) = print(io, typeof(direction));