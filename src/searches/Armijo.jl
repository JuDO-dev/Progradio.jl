struct Armijo{F, I} <: ProgradioSearch{F, I}
    η_A::F          #Armijo slope
    β::F            #backtracking fraction
    k_0::I          #initial backtracking exponent
    k_min::I        #minimum backtracking exponent

    function Armijo(; η_A::F=0.5, β::F=0.6, k_0::I=0, k_min::I=-20) where {F<:AbstractFloat, I<:Integer}
        if !(0 ≤ η_A ≤ 0.5)
            throw(DomainError(η_A, "Ensure that 0 < η_A ≤ ½"))
        end
        if !(zero(F) < β < one(F))
            throw(DomainError(β, "Ensure that 0 < β < 1"))
        end
        if !(k_0 > k_min)
            throw(DomainError(k_0, "Ensure that k_0 > k_min"))
        end
        return new{F, I}(η_A, β, k_0, k_min)
    end
end

mutable struct ArmijoState{F, I} <: ProgradioSearchState{F, I}
    dϕ0::F
    k::I
    α::F
    x_trial::Vector{F}
    α_2::F
    ϕα_2::F
end

function search_state(problem::ProgradioProblem{F, I}, search::Armijo{F, I}) where {F<:AbstractFloat, I<:Integer}
    return ArmijoState(zero(F), search.k_0, zero(F), zeros(F, length(problem.x_0)), zero(F), zero(F))
end

function initiate(problem::UProblem{F, I}, direction::ProgradioDirection{F, I}, search::Armijo{F, I}) where 
    {F<:AbstractFloat, I<:Integer}

    # Initialize state
    n_x = length(problem.x_0);
    state = IteratorState(Iterating(), zero(I),
        deepcopy(problem.x_0), convert(F, Inf), zeros(F, n_x),
        zeros(F, n_x), zero(F), zeros(F, n_x),
        direction_state(problem, direction), search_state(problem, search),
        falses(n_x), falses(n_x), falses(n_x), trues(n_x), #unused
        zero(I), zero(F) #unused
    );
    
    # Evaluate f(x) and ∇f(x)
    state.fx = problem.f(state.x);
    problem.g!(state.gx, state.x);

    return (state.fx, state)
end

function iterate!(state::IteratorState{F, I, DS, SS}, problem::UProblem{F, I}, direction::ProgradioDirection{F, I}, search::Armijo{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Keep previous iteration values in memory
    memorize!(state);

    # Compute a descent direction
    direction!(state, problem, direction);

    # Search for an Armijo step
    state.search_state.dϕ0 = dot(state.gx, state.direction_state.d);
    search!(state, problem, search);

    # Seek improved step via interpolation
    #improve!(state, problem, search);

    # Compute gradient
    problem.g!(state.gx, state.x);

    return nothing
end

function memorize!(state::IteratorState{F, I, DS, SS}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}
    
    @. state.x_previous = state.x;
    state.fx_previous = state.fx;
    return nothing
end

function search!(state::IteratorState{F, I, DS, SS}, problem::ProgradioProblem{F, I}, search::Armijo{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    # Compute initial trial step using previous k
    state.search_state.α, state.fx, admissible = trial_step!(state, problem, search, state.search_state.k, state.x);
    
    # If non-admissible, start back-tracking...
    if !admissible
        # ... until an admissible step is found
        while !admissible
            # Keep previous α and ϕ(α) for interpolation
            state.search_state.α_2 = state.search_state.α;
            state.search_state.ϕα_2 = state.fx;

            # Take trial step
            state.search_state.k += 1;
            state.search_state.α, state.fx, admissible = trial_step!(state, problem, search, state.search_state.k, state.x);
        end

    # Forth-track otherwise
    else
        k = state.search_state.k;
        while true
            # Stop if trial step becomes too large
            if k - 1 < search.k_min
                break
            end

            # Take trial step
            k -= 1;
            α_trial, fx_trial, admissible = trial_step!(state, problem, search, k, state.search_state.x_trial);

            # Stop if trial step becomes non-admissible
            if !admissible
                # Keep trial α and ϕ(α) for interpolation
                state.search_state.α_2 = α_trial;
                state.search_state.ϕα_2 = fx_trial;
                
                break
            end
            
            # Otherwise, keep the trial step
            state.search_state.k = k;
            state.search_state.α = α_trial;
            @. state.x = state.search_state.x_trial;
            state.fx = fx_trial;
        end
    end
    return nothing
end

function trial_step!(state::IteratorState{F, I, DS, SS}, up::UProblem{F, I}, search::Armijo{F, I}, k_trial::I, x_trial::Vector{F}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}
    
    # Take trial step
    α_trial = search.β ^ k_trial;
    @. x_trial = state.x_previous + α_trial * state.direction_state.d;
    fx_trial = up.f(x_trial);

    # Test Armijo condition
    admissible = is_Armijo(state, α_trial, fx_trial, search.η_A);

    return α_trial, fx_trial, admissible
end

function is_Armijo(state::IteratorState{F, I, DS, SS}, α::F, ϕα::F, η_A::F) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}
    
    return ϕα - state.fx_previous ≤ η_A * α * state.search_state.dϕ0
end
#=
function improve!(state::IteratorState{F, I, DS, SS}, problem::ProgradioProblem{F, I}, search::Armijo{F, I}) where 
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
=#
function interpolate(ϕ0::F, α_1::F, ϕ_1::F, α_2::F, ϕ_2::F) where {F<:AbstractFloat}
    a = (α_1 * (ϕ0 - ϕ_2) + α_2 * (ϕ_1 - ϕ_2)) / (α_1 * α_2 * (α_1 - α_2));
    b = (ϕ_1 - ϕ0) / α_1 - a * α_1;
    return -0.5 * b / a
end

# Pretty printing
Base.show(io::IO, search::Armijo) = print(io, typeof(search), " with slope η_A = ", search.η_A, " and fraction β = ", search.β);