# Implementing the CommonSolve interface
function solve(problem::ProgradioProblem{F}, direction::ProgradioDirection{F}, search::ProgradioSearch{F}; 
    i_max::Int=1000, x_tol::F=1e-6, f_tol::F=1e-6, g_tol::F=1e-6) where {F<:AbstractFloat}

    iterator, state = init(problem, direction, search; i_max=i_max, x_tol=x_tol, f_tol=f_tol, g_tol=g_tol);
    return solve!(iterator, state)
end

function init(problem::ProgradioProblem{F}, direction::ProgradioDirection{F}, search::ProgradioSearch{F}; 
    i_max::Int=1000, x_tol::F=1e-6, f_tol::F=1e-6, g_tol::F=1e-6) where {F<:AbstractFloat}

    iterator = Iterator(problem, direction, search; i_max=i_max, x_tol=x_tol, f_tol=f_tol, g_tol=g_tol);
    _, state = iterate(iterator);
    return iterator, state
end

function solve!(iterator::Iterator{F, P, D, S}, state::IteratorState{F, DS, SS}) where 
    {F<:AbstractFloat, P<:ProgradioProblem{F}, D<:ProgradioDirection{F}, S<:ProgradioSearch{F},
    DS<:ProgradioDirectionState{F}, SS<:ProgradioSearchState{F}}

    while true
        next = iterate(iterator, state);
        if next !== nothing
            (_, state) = next;
        else
            break
        end
    end
    return state
end