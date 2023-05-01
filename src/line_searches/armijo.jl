struct Armijo{R<:Real} <: LineSearch
    slope::R
    backtracking::R
    k_min::Int
    k_0::Int
    
    function Armijo(; slope::R=0.5, backtracking::R=0.6, k_0::Integer=0, k_min::Integer=-20) where {R<:Real}

        0 ≤ slope ≤ 0.5         ? nothing : throw(DomainError(slope, "Ensure that 0 < slope ≤ 0.5"));
        0 < backtracking < 1    ? nothing : throw(DomainError(backtracking, "Ensure that 0 < backtracking < 1"));
        k_0 ≥ k_min             ? nothing : throw(DomainError(k_0, "Ensure that k_0 ≥ k_min"));

        return new{R}(slope, backtracking, k_min, k_0)
    end
end

mutable struct ArmijoState{R<:Real, X<:AbstractVector{R}} <: LineSearchState
    dϕ0::R
    k::Int
    α::R
    x_trial::X
    Δx::X
    α_previous::R
    ϕα_previous::R
end

build_state(x::AbstractVector{R}, armijo::Armijo) where {R<:Real} = ArmijoState(zero(R), armijo.k_0,
    zero(R), zero(x), zero(x), zero(R), zero(R)
);

function line_search!(state::IteratorState, problem::Problem, armijo::Armijo)

    state.algorithm_state.line_search_state.dϕ0 = dot(state.gx, state.algorithm_state.direction_state.d, state.set_state.W);
    
    armijo!(state, problem, armijo);

    improve!(state, problem, armijo);

    return nothing
end

function armijo!(state::IteratorState, problem::Problem, armijo::Armijo)

    armijo_state = state.algorithm_state.line_search_state;
    
    armijo_state.α, state.fx, admissible = trial_step!(state, armijo_state.k, state.x, problem, armijo);

    if !admissible
        # Back-track until an admissible step is found
        while !admissible
            armijo_state.α_previous = armijo_state.α;
            armijo_state.ϕα_previous = state.fx;
            armijo_state.k += 1;
            armijo_state.α, state.fx, admissible = trial_step!(state, armijo_state.k, state.x, problem, armijo);
        end

    else
        # Forth-track otherwise
        k = armijo_state.k;
        while true
            # Stop if the trial step grows too large
            if k - 1 < armijo.k_min
                break
            end

            k -= 1;
            α_trial, fx_trial, admissible = trial_step!(state, k, armijo_state.x_trial, problem, armijo);

            # Stop if trial step becomes non-admissible
            if !admissible
                # Keep trial α and ϕ(α) for interpolation
                armijo_state.α_previous = α_trial;
                armijo_state.ϕα_previous = fx_trial;
                
                break
            end
            
            # Otherwise, keep the trial step
            armijo_state.k = k;
            armijo_state.α = α_trial;
            @. state.x = armijo_state.x_trial;
            state.fx = fx_trial;
        end
    end

    return nothing
end

function trial_step!(state::IteratorState, k_trial::Integer, x_trial::AbstractVector, problem::Problem, armijo::Armijo)

    α_trial = armijo.backtracking^k_trial;
    @. x_trial = state.x_previous + α_trial * state.algorithm_state.direction_state.d;
    project!(x_trial, problem.set);
    fx_trial = problem.f(x_trial);

    dϕ0 = state.algorithm_state.line_search_state.dϕ0;
    Δx = state.algorithm_state.line_search_state.Δx;
    @. Δx = state.x_previous - x_trial;
    Δϕ0_B = dot(state.gx, Δx, state.set_state.B);

    admissible = fx_trial - state.fx_previous ≤ armijo.slope * (α_trial * dϕ0 - Δϕ0_B);

    return α_trial, fx_trial, admissible
end

function improve!(state::IteratorState, problem::Problem, armijo::Armijo)

    armijo_state = state.algorithm_state.line_search_state;
    d = state.algorithm_state.direction_state.d;
    
    α_interpolated = interpolate(state.fx_previous, armijo_state.α, state.fx, armijo_state.α_previous, armijo_state.ϕα_previous);

    in_bounds = 0 ≤ α_interpolated ≤ armijo.backtracking ^ armijo.k_min;
    not_redundant = 0 != α_interpolated != armijo_state.α;

    if in_bounds && not_redundant
        if α_interpolated ≥ armijo_state.α
            @. armijo_state.x_trial = state.x_previous + α_interpolated * d;
            
        else
            for j in eachindex(armijo_state.x_trial)
                if state.set_state.W[j]
                    armijo_state.x_trial[j] = state.x_previous[j] + α_interpolated * d[j];
                else
                    armijo_state.x_trial[j] = state.x_previous[j] + armijo_state.α * d[j];
                end
            end
        end

        # Check for improvement in f
        project!(armijo_state.x_trial, problem.set);
        fx_trial = problem.f(armijo_state.x_trial);
        if fx_trial < state.fx
            @. state.x = armijo_state.x_trial;
            state.fx = fx_trial;
        end
    end

    return nothing
end

function interpolate(ϕ0::R, α_1::R, ϕα_1::R, α_2::R, ϕα_2::R) where {R<:Real}

    a = (α_2 * (ϕα_1 - ϕ0) + α_1 * (ϕ0 - ϕα_2));
    b = (α_2 * α_2 * (ϕ0 - ϕα_1) + α_1 * α_1 * (ϕα_2 - ϕ0));
    return -b / (2 * a)
end