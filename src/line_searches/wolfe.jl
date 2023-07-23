struct Wolfe{T<:Real} <: LineSearch
    η_A::T          # Sufficient decrease constant ∈ (0,0.5)
    η_W::T          # Gradient condition constant ∈ (η_A,1)
    γ_e::T          # Step size increase factor ∈ (1,∞)
    α_0::T          # Initial step size guess
    α_max::T        # Max allowable step size

    function Wolfe(; η_A::T=1.0e-4, η_W::T=0.9, γ_e::T=2.0, α_0::T=1.0, α_max::T=5.0) where {T<:Real}

        0 < η_A ≤ 0.5           ? nothing : throw(DomainError(η_A, "Ensure that 0 < η_A ≤ 0.5"));
        η_A < η_W < 1           ? nothing : throw(DomainError(η_W, "Ensure that η_A < η_W < 1"));
        γ_e > 1                 ? nothing : throw(DomainError(γ_e, "Ensure that γ_e > 1"));
        α_0 ≤ α_max             ? nothing : throw(DomainError(α_0, "Ensure that α_0 ≤ α_max"));

        return new{T}(η_A, η_W, γ_e, α_0, α_max)
    end
end

mutable struct WolfeState{T, X<:AbstractArray} <: LineSearchState{T}
    n::Integer
    x_trial::X
    g_trial::X
    d::X
    α::T
    α_old::T
    α_new::T
    κ::X                                # Kink candidate distances in α
    κ_index::Vector{Integer}            # Kink candidate order in α
    κ_range1::UnitRange{Integer}
    κ_range2::UnitRange{Integer}
    ω1::T
    ω2::T
    W1_l::BitVector
    W1_r::BitVector
    dϕα_l::T
    dϕα_r::T
    dϕ0_r::T
end

function build_state(x::AbstractArray{T, N}, wolfe::Wolfe) where {T<:Real, N}
    n = length(x);
    return WolfeState(n, x, deepcopy(x), deepcopy(x), wolfe.α_0, zero(T), zero(T), deepcopy(x), zeros(Integer,n), 
    UnitRange{Integer}(1,n), UnitRange{Integer}(1,n), zero(T), zero(T), BitVector(undef, n), BitVector(undef, n), 
    zero(T), zero(T), zero(T));
end

function line_search!(state::IteratorState, problem::Problem, wolfe::Wolfe)

    state.algorithm_state.line_search_state.dϕ0_r = dot(state.gx, state.algorithm_state.direction_state.d, state.constraints_state.W);
    
    wolfeStageOne!(state, problem, wolfe);

    return nothing
end

function wolfeStageOne!(state::IteratorState, problem::Problem, wolfe::Wolfe) 
    # 1 values correspond to α, 2 values correspond to α_old
    wolfe_state = state.algorithm_state.line_search_state;

    # Obtain line search properties at 0
    ϕ0 = state.fx_previous;
    dϕ0_r = wolfe_state.dϕ0_r;
    wolfe_state.d .= state.algorithm_state.direction_state.d;

    # Reinitialize α
    wolfe_state.α_old = 0.0;

    # ω for α_old is 0 if α_old = 0
    wolfe_state.ω2 = 0.0;

    # Obtain line search properties at α
    wolfe_state.x_trial = state.x_previous + wolfe_state.α * wolfe_state.d;
    project!(wolfe_state.x_trial, problem.constraints);
    ϕα = problem.f(wolfe_state.x_trial);
    wolfe_state.ω1 = ω(wolfe_state.α, ϕ0, ϕα, dϕ0_r, wolfe.η_A);
    problem.g!(wolfe_state.g_trial,wolfe_state.x_trial);

    # Identify kinks and sort them in order of distance
    kinkDist!(wolfe_state.κ, state.x_previous, wolfe_state.d, problem.constraints)
    wolfe_state.κ_index = sortperm(wolfe_state.κ);
    wolfe_state.κ_range1 = searchsorted(wolfe_state.κ[wolfe_state.κ_index],wolfe_state.α);

    # Initialize κ_range2
    wolfe_state.κ_range2 = searchsorted(wolfe_state.κ[wolfe_state.κ_index],wolfe_state.α_old);

    # Compute left and right derivatives at α
    getRLDeriv!(state)

    # Search if α is not a wolfe step and while α is less than α_max
    while !isQWolfe(wolfe_state.ω1, dϕ0_r, wolfe_state.dϕα_l, wolfe_state.dϕα_r, wolfe.η_W) && wolfe_state.α ≠ wolfe.α_max
        if wolfe_state.ω1 ≥ wolfe_state.ω2
            return wolfeStageTwo!(state,problem,wolfe,wolfe_state.α_old,wolfe_state.α,UnitRange{Integer}(clamp(wolfe_state.κ_range2.start,1,wolfe_state.n),clamp(wolfe_state.κ_range1.start,1,wolfe_state.n)),ϕ0,dϕ0_r,false);
        elseif dω(dϕ0_r, wolfe_state.dϕα_l, wolfe.η_A) < 0
            wolfe_state.ω2 = wolfe_state.ω1;
            wolfe_state.κ_range2 .= wolfe_state.κ_range1;
            return wolfeStageTwo!(state,problem,wolfe,wolfe_state.α,wolfe_state.α_old,UnitRange{Integer}(clamp(wolfe_state.κ_range1.start,1,wolfe_state.n),clamp(wolfe_state.κ_range2.start,1,wolfe_state.n)),ϕ0,dϕ0_r,false);
        else
            # Set α_old to α
            state.α_old = state.α;
            state.ω2 = state.ω1;
            state.κ_range2 .= state.κ_range1;

            # Update α
            wolfe_state.α = min(wolfe_state.α*wolfe.γ_e,wolfe.α_max)

            # Obtain line search properties at α
            wolfe_state.x_trial = state.x_previous + wolfe_state.α * wolfe_state.d;
            project!(wolfe_state.x_trial, problem.constraints);
            ϕα = problem.f(wolfe_state.x_trial);
            wolfe_state.ω1 = ω(wolfe_state.α, ϕ0, ϕα, dϕ0_r, wolfe.η_A);
            problem.g!(wolfe_state.g_trial,wolfe_state.x_trial);

            wolfe_state.κ_range1 = searchsorted(wolfe_state.κ[wolfe_state.κ_index],wolfe_state.α);

            # Check if α is a kink point
            getRLDeriv!(state)       
        end
    end

    # Return state with current α if while loop conditions are not satisfied
    @. state.x = wolfe_state.x_trial;
    state.fx = ϕα;
    @. state.gx = wolfe_state.g_trial;
    return nothing
end

function wolfeStageTwo!(state::IteratorState, problem::Problem, wolfe::Wolfe, α_low::T, α_high::T, κ_span::UnitRange{Integer}, ϕ0::T, dϕ0_r::T, callState::Bool) where T<:Real
    # TODO: Add counter for kink reductions before interpolation
    # Note: κ_span covers the range of constraints that switch between α_low and α_high. start is associated with α_low and refers to
    # the index of the constraint (smallest or largest α) nearest to switching at α_low. end serves the same role for α_high.
    
    wolfe_state = state.algorithm_state.line_search_state;

    if callState
        if κ_span.start < κ_span.stop
            # If the interval contains a kink point and α_low is smaller than α_high, set α_new to the next kink point.
            κ_span = UnitRange{Integer}(searchsortedlast(view(wolfe_state.κ[wolfe_state.κ_index],κ_span),α_low) + 1,κ_span.stop);
            wolfe_state.α_new = wolfe_state.κ[wolfe_state.κ_index[κ_span.start]];
        elseif κ_span.start > κ_span.stop
            # If the interval contains a kink point and α_low is larger than α_high, set α_new to the previous kink point.
            κ_span = UnitRange{Integer}(searchsortedfirst(view(wolfe_state.κ[wolfe_state.κ_index],κ_span),α_low) - 1,κ_span.stop);
            wolfe_state.α_new = wolfe_state.κ[wolfe_state.κ_index[κ_span.start]];
        else
            # If the interval contains no kink points, interpolate
            # TODO: Interpolate properly
            wolfe_state.α_new = convert(T,0.5)*(α_high + α_low);
        end
    else
        # Bisect interval on first call or if wolfe counter exceeds limits
        wolfe_state.α_new = convert(T,0.5)*(α_high + α_low);
    end

    # Evaluate left and right derivatives at α_new, as well as the binding set
    wolfe_state.κ_range1 = searchsorted(wolfe_state.κ[wolfe_state.κ_index],wolfe_state.α_new);

    # Obtain line search properties at α_new
    wolfe_state.x_trial = state.x_previous + wolfe_state.α_new * wolfe_state.d;
    project!(wolfe_state.x_trial, problem.constraints);
    ϕα = problem.f(wolfe_state.x_trial);
    wolfe_state.ω1 = ω(wolfe_state.α_new, ϕ0, ϕα, dϕ0_r, wolfe.η_A);
    problem.g!(wolfe_state.g_trial,wolfe_state.x_trial);


    # Check if α is a kink point
    getRLDeriv!(state)

    if !isQWolfe(wolfe_state.ω1, dϕ0_r, wolfe_state.dϕα_l, wolfe_state.dϕα_r, wolfe.η_W)
        if wolfe_state.ω1 ≥ wolfe_state.ω2
            return wolfeStageTwo!(state,problem,wolfe,α_low,wolfe_state.α_new,UnitRange{Integer}(κ_span.start,clamp(wolfe_state.κ_range1.start,1,wolfe_state.n)),ϕ0,dϕ0_r,true);
        elseif wolfe_state.dϕα_r < zero(T) && α_low < α_high
            wolfe_state.ω2 = wolfe_state.ω1;
            wolfe_state.κ_range2 .= wolfe_state.κ_range1;
            return wolfeStageTwo!(state,problem,wolfe,wolfe_state.α_new,α_high,UnitRange{Integer}(clamp(wolfe_state.κ_range1.start,1,wolfe_state.n),κ_span.stop),ϕ0,dϕ0_r,true);
        elseif wolfe_state.dϕα_l > zero(T) && α_low > α_high
            wolfe_state.ω2 = wolfe_state.ω1;
            wolfe_state.κ_range2 .= wolfe_state.κ_range1;
            return wolfeStageTwo!(state,problem,wolfe,wolfe_state.α_new,α_high,UnitRange{Integer}(clamp(wolfe_state.κ_range1.start,1,wolfe_state.n),κ_span.stop),ϕ0,dϕ0_r,true);
        else
            wolfe_state.ω2 = wolfe_state.ω1;
            wolfe_state.κ_range2 .= wolfe_state.κ_range1;
            return wolfeStageTwo!(state,problem,wolfe,wolfe_state.α_new,α_low,UnitRange{Integer}(clamp(wolfe_state.κ_range1.start,1,wolfe_state.n),κ_span.start),ϕ0,dϕ0_r,true);
        end
    else
        wolfe_state.α = wolfe_state.α_new;
        return nothing
    end

end

function ω(α::T, ϕ0::T, ϕα::T, dϕ0_r::T, η_A::T) where T<:Real
    return ϕα - (ϕ0 + α*η_A*dϕ0_r);
end

function dω(dϕ0_r::T, dϕα::T, η_A::T) where T<:Real
    return dϕα - η_A*dϕ0_r;
end

function isQWolfe(ω::T, dϕ0_r::T, dϕα_l::T, dϕα_r::T, η_W::T) where T<:Real
    # Check derivative conditions
    isWolfe1Var = abs(dϕα_l) ≤ η_W * abs(dϕ0_r);
    isWolfe2Var = abs(dϕα_r) ≤ η_W * abs(dϕ0_r);
    isWolfe3Var = (!(dϕα_l ≈ dϕα_r)) && (dϕα_l ≤ 0 ≤ dϕα_r);

    # Return 'true' if quasi-Wolfe conditions are satisfied
    return (ω ≤ 0) && (isWolfe1Var || isWolfe2Var || isWolfe3Var);
end

function getRLDeriv!(state::IteratorState) 
    wolfe_state = state.algorithm_state.line_search_state;
    if wolfe_state.κ_range1.stop < wolfe_state.κ_range1.start
        # Left and right derivatives are equal if α is not a Kink
        binding_from_indices!(wolfe_state.W1_l, view(wolfe_state.κ_index,1:(wolfe_state.κ_range1.start-1)), wolfe_state.n);
        wolfe_state.dϕα_l = dot(wolfe_state.g_trial, wolfe_state.d, wolfe_state.W1_l);
        wolfe_state.W1_r .= wolfe_state.W1_l;
        wolfe_state.dϕα_r = wolfe_state.dϕα_l;
    else
        # Compute left and right derivatives seperately otherwise
        binding_from_indices!(wolfe_state.W1_l, view(wolfe_state.κ_index,1:(wolfe_state.κ_range1.start-1)), wolfe_state.n);
        binding_from_indices!(wolfe_state.W1_r, view(wolfe_state.κ_index,1:wolfe_state.κ_range1.stop), wolfe_state.n);
        wolfe_state.dϕα_l = dot(wolfe_state.g_trial, wolfe_state.d, wolfe_state.W1_l);
        wolfe_state.dϕα_r = dot(wolfe_state.g_trial, wolfe_state.d, wolfe_state.W1_r);
    end
end