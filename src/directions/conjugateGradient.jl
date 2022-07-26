abstract type ConjugateGradient{F} <: ProgradioDirection{F} end 

# State
mutable struct ConjugateGradientState{F} <: ProgradioDirectionState{F}
    r::Integer
    gx_previous::Vector{F}
    d_previous::Vector{F}
end

function direction_state(n::Integer, ::ConjugateGradient{F}) where F<:AbstractFloat
    return ConjugateGradientState{F}(0, ones(F, n), zeros(F, n))
end

# Variants
struct FletcherReeves{F} <: ConjugateGradient{F}
    restart::Integer
    σ_1::F
    σ_2::F
end
FletcherReeves(r::Integer) = FletcherReeves(r, 0.2, 10.0);
FletcherReeves() = FletcherReeves(10);

struct PolakRibiere{F} <: ConjugateGradient{F}
    restart::Integer
    σ_1::F
    σ_2::F
end
PolakRibiere(r::Integer) = PolakRibiere(r, 0.2, 10.0);
PolakRibiere() = PolakRibiere(10);

struct HagerZhang{F} <: ConjugateGradient{F}
    restart::Integer
    σ_1::F
    σ_2::F
end
HagerZhang(r::Integer) = HagerZhang(r, 0.2, 10.0);
HagerZhang() = HagerZhang(10);



# Direction
function direction!(state::ProgradioOptimizerState{F}, direction::ConjugateGradient{F}) where F<:AbstractFloat
    # Save direction
    @. state.d_state.d_previous = state.d;
    
    # Use steepest direction when restarting
    if state.d_state.r >= direction.restart || state.i < 1
        # Use steepest descent direction
        @. state.d = -state.gx;
        state.d_state.r = 0;
    else
        # Use conjugate direction
        conjugate_direction!(state, direction);

        # If unbounded, overwrite to steepest descent
        d_W = view(state.d, state.W);
        gx_W = view(state.gx, state.W);
        toosmall = -dot(d_W, gx_W) <= direction.σ_1 * dot(gx_W, gx_W);
        toolarge = dot(d_W, d_W) >= direction.σ_2 ^ 2 * dot(gx_W, gx_W);
        if toosmall || toolarge 
            @. d_W = -gx_W;
            state.d_state.r = 0;
        else
            state.d_state.r += 1;
        end
    end

    # Save gradient
    @. state.d_state.gx_previous = state.gx;
    return nothing
end

function conjugate_direction!(state::ProgradioOptimizerState{F}, direction::ConjugateGradient{F}) where F<:AbstractFloat
    # Steepest descent for Binding set indices
    gx_B = view(state.gx, state.B);
    d_B = view(state.d, state.B);
    @. d_B = -gx_B;

    # Conjugate direction for Working set indices
    d_W = view(state.d, state.W);
    gx_W = view(state.gx, state.W);
    μ = coefficient!(state, direction);
    d_previous_W = view(state.d_state.d_previous, state.W);
    @. d_W = μ * d_previous_W - gx_W;
    return nothing
end

function coefficient!(state::ProgradioOptimizerState{F}, ::FletcherReeves{F}) where F<:AbstractFloat
    gx_W = view(state.gx, state.W);
    gx_previous_W = view(state.d_state.gx_previous, state.W);
    return dot(gx_W, gx_W) / dot(gx_previous_W, gx_previous_W);
end

function coefficient!(state::ProgradioOptimizerState{F}, ::PolakRibiere{F}) where F<:AbstractFloat
    gx_W = view(state.gx, state.W);
    @. state.x_a = state.gx - state.d_state.gx_previous;
    x_a_W = view(state.x_a, state.W);

    gx_previous_W = view(state.d_state.gx_previous, state.W);
    return dot(gx_W, x_a_W) / dot(gx_previous_W, gx_previous_W)
end

function coefficient!(state::ProgradioOptimizerState{F}, ::HagerZhang{F}) where F<:AbstractFloat
    gx_W = view(state.gx, state.W);
    @. state.x_a = state.gx - state.d_state.gx_previous;
    x_a_W = view(state.x_a, state.W);

    d_previous_W = view(state.d_state.d_previous, state.W);
    deno = dot(d_previous_W, x_a_W);
    
    @. state.x_b = state.x_a - 2 * dot(x_a_W, x_a_W) / deno * state.d_state.d_previous;
    x_b_W = view(state.x_b, state.W);
    return dot(x_b_W, gx_W) / deno
end