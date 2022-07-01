struct SteepestDescent{F} <: ProgradioDirection{F} end
SteepestDescent() = SteepestDescent{Float64}()

function memorize!(::ProgradioIteratorState{F}, ::Nothing) where F<:AbstractFloat
    return nothing
end

function direction!(state::ProgradioIteratorState{F}, ::SteepestDescent{F}) where F<:AbstractFloat
    @. state.d = -state.g;
    return nothing
end