# Implementing julia's iteration interface
struct Iterator{F<:AbstractFloat, I<:Integer, P<:ProgradioProblem{F, I}, D<:ProgradioDirection{F, I}, S<:ProgradioSearch{F, I}}
    problem::P
    direction::D
    search::S
    # Termination
    i_max::I
    x_tol::F
    f_tol::F
    g_tol::F
    
    function Iterator(problem::P, direction::D, search::S; i_max::I=1000, x_tol::F=1e-6, f_tol::F=1e-6, g_tol::F=1e-6) where 
        {F<:AbstractFloat, I<:Integer, P<:ProgradioProblem{F, I}, D<:ProgradioDirection{F, I}, S<:ProgradioSearch{F, I}}
        
        if i_max < zero(I)
            throw(DomainError(i_max, "Maximum number of iterations must not be negative."))
        end
        if x_tol ≤ zero(F) || f_tol ≤ zero(F) || g_tol ≤ zero(F)
            throw(DomainError([x_tol, f_tol, g_tol], "All termination tolerances must be positive."))
        end
        return new{F, I, P, D, S}(problem, direction, search, i_max, x_tol, f_tol, g_tol)
    end
end

Base.IteratorSize(::Iterator) = Base.SizeUnknown();
Base.eltype(::Iterator{F, I, P, D, S}) where {F<:AbstractFloat, I, P, D, S} = F;

struct Iterating end
struct Converged end
struct Exhausted end

mutable struct IteratorState{F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}
    status::Union{Iterating, Converged, Exhausted}
    i::I
    x::Vector{F}
    fx::F
    gx::Vector{F}
    x_previous::Vector{F}
    fx_previous::F
    Δx::Vector{F}
    direction_state::DS
    search_state::SS
    # Box
    B_ℓ::BitVector
    B_u::BitVector
    B::BitVector
    W::BitVector
    # Simplex
    j̄::I
    x_j̄::F
end

# First iteration
function Base.iterate(iterator::Iterator{F, I, P, D, S}) where 
    {F<:AbstractFloat, I<:Integer, P<:ProgradioProblem{F, I}, D<:ProgradioDirection, S<:ProgradioSearch}

    return initiate(iterator.problem, iterator.direction, iterator.search)
end

# Remaining iterations
function Base.iterate(iterator::Iterator{F, I, P, D, S}, state::IteratorState{F, I, DS, SS}) where 
    {F<:AbstractFloat, I<:Integer, P<:ProgradioProblem{F, I}, D<:ProgradioDirection{F, I}, S<:ProgradioSearch{F, I}, 
    DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    if state.i >= iterator.i_max
        state.status = Exhausted();
        return nothing
    elseif state.i > 1 && converged!(state, iterator.x_tol, iterator.f_tol, iterator.g_tol, iterator.problem)
        state.status = Converged();
        return nothing
    else
        iterate!(state, iterator.problem, iterator.direction, iterator.search);
        state.i += one(I);
        return (state.fx, state)
    end
end

# Convergence tests
function converged!(state::IteratorState{F, I, DS, SS}, x_tol::F, f_tol::F, g_tol::F, problem::ProgradioProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    if !f_converged(state, f_tol, problem)
        return false
    elseif !g_converged(state, g_tol, problem)
        return false
    elseif !x_converged!(state, x_tol, problem)
        return false
    else
        return true
    end
end

function x_converged!(state::IteratorState{F, I, DS, SS}, x_tol::F, ::ProgradioProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    @. state.Δx = state.x - state.x_previous;
    norm_Δx = norm(state.Δx);
    return norm_Δx / (1 + norm_Δx) < x_tol ? true : false
end

function f_converged(state::IteratorState{F, I, DS, SS}, f_tol::F, ::ProgradioProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    return abs(state.fx - state.fx_previous) / (1 + abs(state.fx)) < f_tol ? true : false
end

function g_converged(state::IteratorState{F, I, DS, SS}, g_tol::F, ::ProgradioProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    return norm(state.gx, state.W) / (sum(state.W) * (1 + abs(state.fx))) < g_tol ? true : false
end

function g_converged(state::IteratorState{F, I, DS, SS}, g_tol::F, ::UProblem{F, I}) where 
    {F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

    return norm(state.gx) / (length(state.x) * (1 + abs(state.fx))) < g_tol ? true : false
end

# Pretty printing
Base.show(io::IO, ::Iterator{F, I, P, D, S}) where {F, I, P, D, S} = print(io, "Iterator{", F, ",", I, "}");
Base.show(io::IO, state::IteratorState{F, I, DS, SS}) where {F, I, DS, SS} = print(io, "IteratorState{", F, ",", I, "} at iteration ", state.i, " with f(x) = ", state.fx);