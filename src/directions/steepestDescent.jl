struct SteepestDescent{F} <: ProgradioDirection{F} end
SteepestDescent() = SteepestDescent{Float64}()

struct SteepestDescentState{F} <: ProgradioDirectionState{F} end
direction_state(::Integer, ::SteepestDescent{F}) where F = SteepestDescentState{F}();

function direction!(state::ProgradioOptimizerState{F}, ::SteepestDescent{F}) where F<:AbstractFloat
    @. state.d = -state.gx;
    return nothing
end