function initiate(problem::BCProblem{F, I}, direction::ProgradioDirection{F, I}, search::Armijo{F, I}) where 
    {F<:AbstractFloat, I<:Integer}

    # Initialize state
    n_x = length(problem.x_0);
    state = IteratorState(Iterating(), zero(I),
        deepcopy(problem.x_0), convert(F, Inf), zeros(F, n_x),
        zeros(F, n_x), zero(F), zeros(F, n_x),
        direction_state(problem, direction), search_state(problem, search),
        falses(n_x), falses(n_x), falses(n_x), trues(n_x),
        zero(I), zero(F) #unused
    );
    
    # Evaluate f(x_0) and ∇f(x_0)
    state.fx = problem.f(state.x);
    problem.g!(state.gx, state.x);

    # Compute binding sets
    binding!(state, problem);

    return (state.fx, state)
end

function binding!(state::IteratorState{F, I, DS, SS}, problem::BCProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState, SS<:ProgradioSearchState}

    # Binding tolerance
    @. state.Δx = state.x - state.gx;
    project!(state.Δx, problem.ℓ, problem.u);
    @. state.Δx = state.x - state.Δx;
    ϵ = min(norm(state.Δx), problem.B_tol);

    # Identify binding indices
    for j in eachindex(state.x)
        if problem.ℓ[j] ≤ state.x[j] ≤ (problem.ℓ[j] + ϵ) && state.gx[j] > 0
            state.B_ℓ[j] = true;
            state.B_u[j] = false;
            state.B[j] = true;
            state.W[j] = false;
        elseif (problem.u[j] - ϵ) ≤ state.x[j] ≤ problem.u[j]  && state.gx[j] < 0
            state.B_ℓ[j] = false;
            state.B_u[j] = true;
            state.B[j] = true;
            state.W[j] = false;
        else
            state.B_ℓ[j] = false;
            state.B_u[j] = false;
            state.B[j] = false;
            state.W[j] = true;
        end
    end
    return nothing
end

function iterate!(state::IteratorState{F, I, DS, SS}, problem::BCProblem{F, I}, direction::ProgradioDirection{F, I}, search::Armijo{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Keep previous iteration values in memory
    memorize!(state);

    # Compute a descent direction
    direction!(state, problem, direction);
    
    # Projected search for an Armijo-like step
    state.search_state.dϕ0 = dot(state.gx, state.direction_state.d, state.W);
    search!(state, problem, search);
    
    # Seek improved step via interpolation
    improve!(state, problem, search);
    
    # Compute gradient
    problem.g!(state.gx, state.x);
        
    # Update binding sets
    binding!(state, problem);
    
    return nothing
end

function trial_step!(state::IteratorState{F, I, DS, SS}, problem::BCProblem{F, I}, search::Armijo{F, I}, k_trial::I, x_trial::Vector{F}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Take trial step 
    α_trial = search.β ^ k_trial;
    @. x_trial = state.x_previous + α_trial * state.direction_state.d;
    project!(x_trial, problem.ℓ, problem.u);
    fx_trial = problem.f(x_trial);

    # Δϕ(0) for binding set
    @. state.Δx = state.x_previous - x_trial;
    Δϕ0_B = dot(state.gx, state.Δx, state.B);

    # Test Armijo condition
    admissible = is_Armijo(state, α_trial, fx_trial, search.η_A, Δϕ0_B);

    return α_trial, fx_trial, admissible
end

function is_Armijo(state::IteratorState{F, I, DS, SS}, α::F, ϕα::F, η_A::F, Δϕ0_B::F) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}
    
    return ϕα - state.fx_previous ≤ η_A * (α * state.search_state.dϕ0 - Δϕ0_B)
end

function improve!(state::IteratorState{F, I, DS, SS}, problem::BCProblem{F, I}, search::Armijo{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Interpolate a new α using the two previous line-search iterations
    α_trial = interpolate(state.fx_previous, state.search_state.α, state.fx, state.search_state.α_2, state.search_state.ϕα_2);

    # Check if bounded
    if 0 < α_trial < search.β ^ search.k_min
        
        # Compute x_trial
        if α_trial ≥ state.search_state.α
            @. state.search_state.x_trial = state.x_previous + α_trial * state.direction_state.d;
        else
            for j in eachindex(state.search_state.x_trial)
                if state.W[j]
                    state.search_state.x_trial[j] = state.x_previous[j] + α_trial * state.direction_state.d[j];
                else
                    state.search_state.x_trial[j] = state.x_previous[j] + state.search_state.α * state.direction_state.d[j];
                end
            end
        end

        # Project and check for improvement in f
        project!(state.search_state.x_trial, problem.ℓ, problem.u);
        fx_trial = problem.f(state.search_state.x_trial);
        if fx_trial < state.fx
            @. state.x = state.search_state.x_trial;
            state.fx = fx_trial;
        end
    end

    return nothing
end