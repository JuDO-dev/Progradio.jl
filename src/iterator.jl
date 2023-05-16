struct Iterator{T<:Real, P<:Problem, A<:Algorithm}
    problem::P
    algorithm::A
    i_max::Int
    g_tol::T

    function Iterator(problem::P, algorithm::A; i_max::Integer=10_000, g_tol::T=1e-8) where {T<:Real,
        X, C, P<:Problem{T, X, C}, A<:Algorithm}

        i_max > 0 ? nothing : throw(DomainError(i_max, "Ensure i_max > 0"));
        g_tol > 0 ? nothing : throw(DomainError(g_tol, "Ensure g_tol > 0"));
        
        return new{T, P, A}(problem, algorithm, i_max, g_tol)
    end
end

Base.eltype(::Iterator{T, P, A}) where {T, P, A} = T;
Base.IteratorSize(::Iterator) = Base.SizeUnknown();

@enum Status iterating converged exhausted;

mutable struct IteratorState{T<:Real, X<:AbstractArray, SS<:ConstraintsState, AS<:AlgorithmState}
    status::Status
    i::Int
    x::X
    fx::T
    gx::Vector{T}
    constraints_state::SS
    algorithm_state::AS
    x_previous::X
    fx_previous::T
end

function build_state(iterator::Iterator{T, P, A}) where {T, P, A}

    x_start = iterator.problem.x_start;
    if !is_inside(x_start, iterator.problem.constraints)
        project!(x_start, iterator.problem.constraints);
    end

    return IteratorState(iterating, 0, deepcopy(x_start), typemax(T), zeros(length(x_start)),
        build_state(x_start, iterator.problem.constraints), build_state(x_start, iterator.algorithm),
        deepcopy(x_start), typemax(T)
    );
end

function Base.iterate(iterator::Iterator)

    state = build_state(iterator);
    
    state.fx = iterator.problem.f(state.x);
    iterator.problem.g!(state.gx, state.x);
    binding!(state, iterator.problem.constraints);

    return (state.fx, state) 
end

function Base.iterate(iterator::Iterator, state::IteratorState)

    if state.i >= iterator.i_max
        state.status = exhausted;
        return nothing

    elseif has_converged(state, iterator.g_tol, iterator.problem.constraints)
        state.status = converged;
        return nothing

    else 
        memorize!(state);
        iterate_algorithm!(state, iterator.problem, iterator.algorithm);
        state.i += 1;
        return (state.fx, state)
    end
end

function memorize!(state::IteratorState)
    
    @. state.x_previous = state.x;
    state.fx_previous = state.fx;
    return nothing
end