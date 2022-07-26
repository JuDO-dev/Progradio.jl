# Optimizer
struct Armijo{F, D} <: ProgradioOptimizer{F, D}
    direction::D
    η_A::F
    β::F
    k_min::Integer
    k_0::Integer
    α_ratio_min::F
end

function Armijo(direction::ProgradioDirection{F}) where F<:AbstractFloat
    return Armijo(direction, convert(F, 0.5), convert(F, 0.6), -20, 0, convert(F, 0.2))
end

# Iterator state
mutable struct ArmijoState{F, DS<:ProgradioDirectionState} <: ProgradioOptimizerState{F, DS}
    i::Integer
    # Memory
    x_previous::Vector{F}
    fx_previous::F
    # Direction
    d::Vector{F}
    d_state::DS
    # Search
    k::Integer
    α::F
    x::Vector{F}
    fx::F
    # Binding
    gx::Vector{F}
    B::BitVector
    W::BitVector
    # Storage
    x_a::Vector{F}
    x_b::Vector{F}
    x_c::Vector{F}
end

function optimizer_state(bci::BCIterator{F, O}) where {F<:AbstractFloat, O<:Armijo}
    x_0 = bci.bcp.x_0;
    n = length(x_0);
    return ArmijoState(0, similar(x_0), zero(F),
        zeros(F, n), direction_state(n, bci.optimizer.direction),
        bci.optimizer.k_0, zero(F), similar(x_0), one(F),
        similar(x_0), BitVector(undef, n), BitVector(undef, n),
        similar(x_0), similar(x_0), similar(x_0))
end

function iterate!(bci::BCIterator, state::ArmijoState)
    # Keep previous iteration data in memory
    memorize!(state);
    
    # Compute a descent direction
    direction!(state, bci.optimizer.direction);

    # Search along this direction
    search!(state, bci.bcp, bci.optimizer);

    # Binding set
    binding!(state, bci.bcp);

    state.i += 1;
    return nothing
end

function memorize!(state::ArmijoState{F}) where F<:AbstractFloat
    @. state.x_previous = state.x;
    state.fx_previous = state.fx;
    return nothing
end

function search!(state::ArmijoState, bcp::BCProblem, optimizer::Armijo)
    # Search until Armijo condition is "just" satisfied
    backforthtrack!(state, bcp, optimizer);

    # Try to improve the step
    improve!(state, bcp, optimizer);

    return nothing
end

function backforthtrack!(state::ArmijoState{F}, bcp::BCProblem{F}, optimizer::Armijo) where F<:AbstractFloat
    # Constant values
    ϕ0 = state.fx_previous;
    dϕ0_W = dot(view(state.gx, state.W), view(state.d, state.W));

    # Start with previous k
    k = state.k;
    α = optimizer.β ^ k;
    @. state.x = projection(state.x_previous + α * state.d, bcp.x_ℓ, bcp.x_u);
    ϕα = bcp.f(state.x);
    
    @. state.x_c = state.x_previous - state.x;
    dϕ0_B = dot(view(state.gx, state.B), view(state.x_c, state.B));

    # Start back-tracking if Armijo condition is not satisfied
    if !isArmijo(α, ϕα, ϕ0, optimizer.η_A, dϕ0_W, dϕ0_B);
        while true
            k += 1;
            α = optimizer.β ^ k;
            @. state.x = projection(state.x_previous + α * state.d, bcp.x_ℓ, bcp.x_u);
            ϕα = bcp.f(state.x);

            @. state.x_c = state.x_previous - state.x;
            dϕ0_B = dot(view(state.gx, state.B), view(state.x_c, state.B));

            # Stop back-tracking
            if isArmijo(α, ϕα, ϕ0, optimizer.η_A, dϕ0_W, dϕ0_B)
                state.k = k;
                state.α = α;
                state.fx = ϕα;
                break
            end
        end
    # Start forth-tracking otherwise
    else
        while true
            # Enforce maximum step-size
            if k <= optimizer.k_min
                break
            end

            k -= 1;
            α_b = optimizer.β ^ k;
            @. state.x_b = projection(state.x_previous + α_b * state.d, bcp.x_ℓ, bcp.x_u);
            ϕα_b = bcp.f(state.x_b);

            @. state.x_c = state.x_previous - state.x_b;
            dϕ0_B = dot(view(state.gx, state.B), view(state.x_c, state.B));

            # Stop forth-tracking
            if !isArmijo(α_b, ϕα_b, ϕ0, optimizer.η_A, dϕ0_W, dϕ0_B)
                break
            else
                state.k = k;
                state.α = α_b;
                @. state.x = state.x_b;
                state.fx = ϕα_b;
            end
        end
    end
    return nothing
end

function isArmijo(α::F, ϕα::F, ϕ0::F, η_A::F, dϕ0_W::F, dϕ0_B::F) where F<:AbstractFloat
    return ϕα - ϕ0 <= η_A * (α * dϕ0_W - dϕ0_B)
end

function improve!(state::ArmijoState{F}, bcp::BCProblem{F}, optimizer::Armijo) where F<:AbstractFloat
    # Seek improved step via interpolation
    @. state.x_b = state.x - state.x_previous;
    dϕ0 = dot(state.gx, state.x_b) / norm(state.x_b);
    α_b = interpolate(state.fx_previous, dϕ0, state.α, state.fx);

    # If interpolated step is within bounds
    if (0 <= α_b <= optimizer.β ^ optimizer.k_min)
        # Check for improvement in f
        α_min = optimizer.α_ratio_min * state.α;
        if α_b >= α_min
            @. state.x_b = projection(state.x_previous + α_b * state.d, bcp.x_ℓ, bcp.x_u);
        else
            x_b_B = view(state.x_b, state.B);
            @. x_b_B = projection(view(state.x_previous, state.B) + α_min * view(state.d, state.B), view(bcp.ℓ, state.B), view(bcp.u, state.B));
            x_b_W = view(state.x_b, state.W);
            @. x_b_W = projection(view(state.x_previous, state.W) + α_b * view(state.d, state.W), view(bcp.ℓ, state.W), view(bcp.u, state.W));
        end
        ϕα_b = bcp.f(state.x_b);
        if ϕα_b <= state.fx
            # Use interpolated step
            @. state.x = state.x_b;
            state.fx = ϕα_b;
        end
    end
    
    return nothing
end

function interpolate(ϕ0::F, dϕ0::F, α::F, ϕα::F) where F<:AbstractFloat
    return -(dϕ0 * α * α) / (2 * (ϕα - ϕ0 - dϕ0 * α))
end