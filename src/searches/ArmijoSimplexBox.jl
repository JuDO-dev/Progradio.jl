function initiate(problem::SBCProblem{F, I}, direction::ProgradioDirection{F, I}, search::Armijo{F, I}) where 
    {F<:AbstractFloat, I<:Integer}

    # Initialize state
    n_x = length(problem.x_0);
    state = IteratorState(Iterating(), zero(I),
        deepcopy(problem.x_0), convert(F, Inf), zeros(F, n_x),
        zeros(F, n_x), zero(F), zeros(F, n_x),
        direction_state(problem, direction), search_state(problem, search),
        falses(n_x), falses(n_x), falses(n_x), trues(n_x),
        zero(I), zero(F)
    );

    # Evaluate f(x_0) and ∇f(x_0)
    state.fx = problem.f(state.x);
    problem.g!(state.gx, state.x);

    # Local coordinate transformation
    remove_simplex!(state, problem);

    # Compute binding sets
    binding!(state, problem);

    # Revert coordinate transformation for x
    state.x[state.j̄] = state.x_j̄;

    return (state.fx, state)
end

function remove_simplex!(state::IteratorState{F, I, DS, SS}, problem::SBCProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Find largest simplex component
    state.x_j̄ = zero(F);
    for j in problem.S
        if state.x[j] > state.x_j̄
            state.j̄ = j;
            state.x_j̄ = state.x[j];
        end
    end
    
    # Transform x
    state.x[state.j̄] = one(F);

    # Transform gx
    for j in problem.S
        if j != state.j̄
            state.gx[j] -= state.gx[state.j̄];
        end
    end

    return nothing
end

function binding!(state::IteratorState{F, I, DS, SS}, problem::SBCProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Update simplex and box binding tolerances
    @. state.Δx = state.x - state.gx;
    project!(state.Δx, problem.S);
    project!(state.Δx, problem.ℓ, problem.u, problem.S̄)
    @. state.Δx = state.x - state.Δx;
    state.Δx[state.j̄] = zero(F);
    S_tol = min(norm(state.Δx, problem.S), problem.S_tol, convert(F, 0.1));
    B_tol = min(norm(state.Δx, problem.S̄), problem.B_tol);

    # Update binding sets
    for j in problem.S
        if j == state.j̄ #x_j̄ bound to 1.0
            state.B_ℓ[state.j̄] = true;
            state.B_u[state.j̄] = true;
            state.B[state.j̄] = true;
            state.W[state.j̄] = false;
        elseif (0 ≤ state.x[j] ≤ S_tol) && (state.gx[j] > 0)
            state.B_ℓ[j] = true;
            state.B_u[j] = false;
            state.B[j] = true;
            state.W[j] = false;
        else
            state.B_ℓ[j] = false;
            state.B_u[j] = false;
            state.B[j] = false;
            state.W[j] = true;
        end
    end
    for j in problem.S̄
        if problem.ℓ[j] ≤ state.x[j] ≤ (problem.ℓ[j] + B_tol) && state.gx[j] > 0
            state.B_ℓ[j] = true;
            state.B_u[j] = false;
            state.B[j] = true;
            state.W[j] = false;
        elseif (problem.u[j] - B_tol) ≤ state.x[j] ≤ problem.u[j]  && state.gx[j] < 0
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


function iterate!(state::IteratorState{F, I, DS, SS}, problem::SBCProblem{F, I}, direction::ProgradioDirection{F, I}, search::Armijo{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Local coordinate transformation for x
    state.x[state.j̄] = one(F);
    
    # Keep previous iteration values in memory
    memorize!(state);

    # Compute a descent direction
    direction!(state, direction);

    # Projected search for an Armijo-like step
    state.search_state.dϕ0 = dot(state.gx, state.direction_state.d, state.W);
    search!(state, problem, search);
    
    # Seek improved step via interpolation
    #improve!(state, bcp, optimizer);

    # Revert coordinate transformation
    revert_simplex!(state, problem);

    # Compute gradient
    problem.g!(state.gx, state.x);

    # Local coordinate transformation
    remove_simplex!(state, problem);

    # Compute binding sets
    binding!(state, problem);

    # Revert coordinate transformation for x
    state.x[state.j̄] = state.x_j̄;

    return nothing
end

function trial_step!(state::IteratorState{F, I, DS, SS}, problem::SBCProblem{F, I}, search::Armijo{F, I}, k_trial::I, x_trial::Vector{F}) where
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Take trial step
    α_trial = search.β ^ k_trial;
    @. x_trial = state.x_previous + α_trial * state.direction_state.d;
    project!(x_trial, problem.S);
    project!(x_trial, problem.ℓ, problem.u, problem.S̄);
    x_trial[state.j̄] = one(F);
    
    # Test simplex membership
    sum_simplex = zero(F);
    for j in problem.S
        if j != state.j̄
            sum_simplex += x_trial[j];
        end
    end

    # Continue if inside simplex
    if sum_simplex > 1
        return α_trial, zero(F), false
    else
        # Evaluate objective [WIP]
        x_trial[state.j̄] -= sum_simplex;
        fx_trial = problem.f(x_trial);
        x_trial[state.j̄] = one(F);

        # Δϕ(0) for binding set
        @. state.Δx = state.x_previous - x_trial;
        Δϕ0_B = dot(state.gx, state.Δx, state.B);

        # Test Armijo condition
        admissible = is_Armijo(state, α_trial, fx_trial, search.η_A, Δϕ0_B);
        
        return α_trial, fx_trial, admissible
    end
end

function revert_simplex!(state::IteratorState{F, I, DS, SS}, problem::SBCProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Transform x and x_previous
    sum_x = zero(F);
    sum_x_previous = zero(F);
    for j in problem.S
        if j != state.j̄
            sum_x += state.x[j];
            sum_x_previous += state.x_previous[j];
        end
    end
    state.x_j̄ = state.x[state.j̄] - sum_x;
    state.x[state.j̄] = state.x_j̄;
    state.x_previous[state.j̄] -= sum_x_previous;
    return nothing
end