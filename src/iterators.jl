abstract type ProgradioIterator{F<:AbstractFloat} end
Base.length(I::ProgradioIterator) = 1 + I.i_max;
Base.eltype(::ProgradioIterator{F}) where F = F;

# Box-constrained problems
struct BCIterator{F, O<:ProgradioOptimizer} <: ProgradioIterator{F}
    # Main
    bcp::BCProblem
    optimizer::O
    # Termination
    i_max::Integer
    f_tol::F
    g_tol::F
    x_tol::F
end

function iterator(bcp::BCProblem{F}, optimizer::O; i_max=20, f_tol=1e-6, g_tol=1e-6, x_tol=1e-6) where {F<:AbstractFloat, O<:ProgradioOptimizer}
    return BCIterator{F, O}(bcp, optimizer, i_max, f_tol, g_tol, x_tol)
end

function Base.iterate(bci::BCIterator)
    # Empty state
    state = optimizer_state(bci);

    # Initial guess
    @. state.x = projection(bci.bcp.x_0, bci.bcp.x_ℓ, bci.bcp.x_u);
    state.fx = bci.bcp.f(state.x);
    
    # Binding set
    binding!(state, bci.bcp);

    return (state.fx, state)
end

function Base.iterate(bci::BCIterator, state::ProgradioOptimizerState{F, DS}) where {F<:AbstractFloat, DS}
    if state.i >= bci.i_max
        return nothing
    elseif state.i > 1 && converged!(state, bci.f_tol, bci.g_tol, bci.x_tol);
        return nothing
    else
        iterate!(bci, state);
        return (state.fx, state)
    end
end

# Sub-routines
function binding!(state::ProgradioOptimizerState{F, DS}, bcp::BCProblem{F}) where {F<:AbstractFloat, DS}
    # Gradient
    bcp.g!(state.gx, state.x);
    # Binding tolerance
    @. state.x_a = state.x - projection(state.x - state.gx, bcp.x_ℓ, bcp.x_u);
    ϵ = min(norm(state.x_a), bcp.ϵ);
    # Set identification
    @. state.B = binding(state.x, bcp.x_ℓ, bcp.x_u, ϵ, state.gx);
    @. state.W = !state.B;
    return nothing
end

function converged!(state::ProgradioOptimizerState{F, DS}, f_tol::F, g_tol::F, x_tol::F) where {F<:AbstractFloat, DS}
    # Objective
    if abs(state.fx - state.fx_previous) / (1 + abs(state.fx)) >= f_tol
        return false
    # Gradient
    elseif norm(view(state.gx, state.W)) / (sum(state.W) * (1 + abs(state.fx))) >= g_tol
        return false
    # Variables
    else
        @. state.x_a = state.x - state.x_previous;
        if norm(state.x_a, Inf) / (1 + norm(state.x, Inf)) >= x_tol
            return false
        else
            return true
        end
    end
end