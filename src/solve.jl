# Implementing the CommonSolve interface
function solve(problem::ProgradioProblem{F, I}, direction::ProgradioDirection{F, I}, search::ProgradioSearch{F, I}; 
    i_max::I=1000, x_tol::F=1e-6, f_tol::F=1e-6, g_tol::F=1e-6) where {F<:AbstractFloat, I<:Integer}

    iterator, state = init(problem, direction, search; i_max=i_max, x_tol=x_tol, f_tol=f_tol, g_tol=g_tol);
    return solve!(iterator, state)
end

function init(problem::ProgradioProblem{F, I}, direction::ProgradioDirection{F, I}, search::ProgradioSearch{F, I}; 
    i_max::I=1000, x_tol::F=1e-6, f_tol::F=1e-6, g_tol::F=1e-6) where {F<:AbstractFloat, I<:Integer}

    iterator = Iterator(problem, direction, search; i_max=i_max, x_tol=x_tol, f_tol=f_tol, g_tol=g_tol);
    _, state = iterate(iterator);
    return iterator, state
end

function solve!(iterator::Iterator{F, I, P, D, S}, state::IteratorState{F, I, DS, SS}) where 
    {F<:AbstractFloat, I<:Integer, P<:ProgradioProblem{F, I}, D<:ProgradioDirection{F, I}, S<:ProgradioSearch{F, I},
    DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}}

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