abstract type CGVariant end
struct FletcherReeves <: CGVariant end
struct PolakRibiere <: CGVariant end
struct HagerZhang <: CGVariant end

struct ConjugateGradient{F, V} <: ProgradioDirection{F} where {V<:CGVariant}
    restart::Int
    σ_min::F
    σ_max::F
    variant::V

    function ConjugateGradient(restart::Int, σ_min::F, σ_max::F, variant::V) where
        {F<:AbstractFloat, V<:CGVariant}

        # Enforce domains
        (restart > 1) ? nothing : throw(DomainError(restart, "Ensure restart > 1"));
        (0 < σ_min < 1) ? nothing : throw(DomainError(σ_min, "Ensure 0 < σ_mim < 1"));
        (σ_max > 1) ? nothing : throw(DomainError(σ_max, "Ensure σ_max > 1"));

        return new{F, V}(restart, σ_min, σ_max, variant)
    end
end

mutable struct ConjugateGradientState{F} <: ProgradioDirectionState{F}
    d::Vector{F}
    d_previous::Vector{F}
    gx_previous::Vector{F}
    r::Int
end

function direction_state(problem::ProgradioProblem{F}, ::ConjugateGradient{F, V}) where
    {F<:AbstractFloat, V<:CGVariant}

    n_x = length(problem.x_0);
    return ConjugateGradientState(zeros(F, n_x), zeros(F, n_x), zeros(F, n_x), 0)
end

function direction!(state::IteratorState{F, DS, SS}, ::UProblem{F}, direction::ConjugateGradient{F, V}) where
    {F<:AbstractFloat, V<:CGVariant, DS<:ConjugateGradientState{F}, SS<:ProgradioSearchState{F}}

    # Keep previous direction in memory
    @. state.direction_state.d_previous = state.direction_state.d;

    # Use steepest descent when restarting
    if state.direction_state.r ≥ direction.restart || state.i < 1
        @. state.direction_state.d = -state.gx;
        state.direction_state.r = 0;
    else
        μ = coefficient!(state, direction.variant);
        μ = isnan(μ) || isinf(μ) ? 0.0 : μ;

        # Compute conjugate direction
        @. state.direction_state.d = μ * state.direction_state.d_previous - state.gx;

        # Revert to steepest descent if conjugate direction is unbounded
        gx_dotted = dot(state.gx, state.gx);
        too_small = -dot(state.direction_state.d, state.gx) ≤ direction.σ_min * gx_dotted;
        too_large = dot(state.direction_state.d, state.direction_state.d) ≥ direction.σ_max * gx_dotted;
        if too_small || too_large    
            @. state.direction_state.d = -state.gx
            state.direction_state.r = 0;
        else
            state.direction_state.r += 1;
        end
    end

    # Keep previous gradient in memory
    @. state.direction_state.gx_previous = state.gx;

    return nothing
end

function direction!(state::IteratorState{F, DS, SS}, ::BCProblem{F}, direction::ConjugateGradient{F, V}) where
    {F<:AbstractFloat, V<:CGVariant, DS<:ConjugateGradientState{F}, SS<:ProgradioSearchState{F}}

    # Keep previous direction in memory
    @. state.direction_state.d_previous = state.direction_state.d;

    # Use steepest descent when restarting
    if state.direction_state.r ≥ direction.restart || state.i < 1
        @. state.direction_state.d = -state.gx;
        state.direction_state.r = 0;
    else
        μ = coefficient!(state, direction.variant);
        μ = isnan(μ) || isinf(μ) ? 0.0 : μ;

        # Compute conjugate direction for W indices
        for j in eachindex(state.direction_state.d, state.W)
            if state.W[j]
                state.direction_state.d[j] = μ * state.direction_state.d_previous[j] - state.gx[j];
            else
                state.direction_state.d[j] = -state.gx[j];
            end
        end
            
        # Revert to steepest descent if conjugate direction is unbounded
        gx_dotted_W = dot(state.gx, state.gx, state.W);
        too_small = -dot(state.direction_state.d, state.gx, state.W) ≤ direction.σ_min * gx_dotted_W;
        too_large = dot(state.direction_state.d, state.direction_state.d, state.W) ≥ direction.σ_max * gx_dotted_W;
        if too_small || too_large    
            for j in eachindex(state.direction_state.d, state.W)
                if state.W[j]
                    state.direction_state.d[j] = -state.gx[j];
                end
            end
            state.direction_state.r = 0;
        else
            state.direction_state.r += 1;
        end
    end

    # Keep previous gradient in memory
    @. state.direction_state.gx_previous = state.gx;

    return nothing
end

function direction!(state::IteratorState{F, DS, SS}, problem::SBCProblem{F}, direction::ConjugateGradient{F, V}) where
    {F<:AbstractFloat, V<:CGVariant, DS<:ConjugateGradientState{F}, SS<:ProgradioSearchState{F}}

    # Keep previous direction in memory
    @. state.direction_state.d_previous = state.direction_state.d;

    # Use steepest descent when restarting
    if state.direction_state.r ≥ direction.restart || state.i < 1
        @. state.direction_state.d = -state.gx;
        state.direction_state.r = 0;
    else
        μ = coefficient!(state, direction.variant);
        μ = isnan(μ) || isinf(μ) ? 0.0 : μ;

        # Compute conjugate direction for W indices of S̄
        for j in problem.S̄
            if state.W[j]
                state.direction_state.d[j] = μ * state.direction_state.d_previous[j] - state.gx[j];
            else
                state.direction_state.d[j] = -state.gx[j];
            end
        end
            
        # Revert to steepest descent if conjugate direction is unbounded
        gx_dotted_W = dot(state.gx, state.gx, state.W);
        too_small = -dot(state.direction_state.d, state.gx, state.W) ≤ direction.σ_min * gx_dotted_W;
        too_large = dot(state.direction_state.d, state.direction_state.d, state.W) ≥ direction.σ_max * gx_dotted_W;
        if too_small || too_large    
            for j in eachindex(state.direction_state.d, state.W)
                if state.W[j]
                    state.direction_state.d[j] = -state.gx[j];
                end
            end
            state.direction_state.r = 0;
        else
            state.direction_state.r += 1;
        end

        # Use steepest descent for S indices
        for j in problem.S
            state.direction_state.d[j] = -state.gx[j];
        end
    end

    # Keep previous gradient in memory
    @. state.direction_state.gx_previous = state.gx;

    return nothing
end

function coefficient!(state::IteratorState{F, DS, SS}, ::FletcherReeves) where
    {F<:AbstractFloat, DS<:ConjugateGradientState{F}, SS<:ProgradioSearchState{F}}

    return dot(state.gx, state.gx, state.W) / dot(state.direction_state.gx_previous, state.direction_state.gx_previous, state.W);
end
    
function coefficient!(state::IteratorState{F, DS, SS}, ::PolakRibiere) where
    {F<:AbstractFloat, DS<:ConjugateGradientState{F}, SS<:ProgradioSearchState{F}}
    
    @. state.Δx = state.gx - state.direction_state.gx_previous;
    return dot(state.gx, state.Δx, state.W) / dot(state.direction_state.gx_previous, state.direction_state.gx_previous, state.W);
end

function coefficient!(state::IteratorState{F, DS, SS}, ::HagerZhang) where
    {F<:AbstractFloat, DS<:ConjugateGradientState{F}, SS<:ProgradioSearchState{F}}

    @. state.Δx = state.gx - state.direction_state.gx_previous;
    denominator = dot(state.direction_state.d_previous, state.Δx, state.W);
    fraction = dot(state.Δx, state.Δx, state.W) / denominator;
    @. state.Δx -= 2 * state.direction_state.d_previous * fraction;
    return dot(state.Δx, state.gx, state.W) / denominator 
end

# Pretty printing
Base.show(io::IO, direction::ConjugateGradient) = print(io, typeof(direction), " with restarts every ", direction.restart, " iterations");