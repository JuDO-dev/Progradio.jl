struct Iterator{R<:Real, P<:Problem, A<:Algorithm}
    problem::P
    algorithm::A
    i_max::Int
    g_tol::R

    function Iterator(problem::P, algorithm::A; i_max::Integer=10_000, g_tol::R=1e-8) where {R, P, A}

        i_max > 0 ? nothing : throw(DomainError(i_max, "Ensure"));
        g_tol > 0 ? nothing : throw(DomainError(g_tol, "Ensure"));
        
        return new{R, P, A}(problem, algorithm, i_max, g_tol)
    end
end

Base.eltype(::Iterator{R, P, A}) where {R, P, A} = R;
Base.IteratorSize(::Iterator) = Base.SizeUnknown();

@enum Status iterating converged exhausted;

mutable struct IteratorState{R<:Real, X<:AbstractVector{R}, SS<:ConstraintSetState, AS<:AlgorithmState}
    status::Status
    i::Int
    x::X
    fx::R
    gx::X
    set_state::SS
    algorithm_state::AS
    x_previous::X
    fx_previous::R
end

function build_state(iterator::Iterator{R, P, A}) where {R, P, A}

    x_guess = iterator.problem.x_guess;

    return IteratorState(iterating, 0, deepcopy(x_guess), typemax(R), zero(x_guess),
        build_state(x_guess, iterator.problem.set), build_state(x_guess, iterator.algorithm),
        zero(x_guess), typemax(R)
    );
end

function Base.iterate(iterator::Iterator)

    state = build_state(iterator);
    
    state.fx = iterator.problem.f(state.x);
    iterator.problem.g!(state.gx, state.x);
    binding!(state, iterator.problem.set);

    return (state.fx, state) 
end

function Base.iterate(iterator::Iterator, state::IteratorState)

    if state.i >= iterator.i_max
        state.status = exhausted;
        return nothing

    elseif has_converged(state, iterator.g_tol, iterator.problem.set)
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