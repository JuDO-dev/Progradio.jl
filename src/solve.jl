function init(problem::Problem, algorithm::Algorithm; i_max::Integer=10_000, g_tol::Real=1e-8)

    iterator = Iterator(problem, algorithm; i_max, g_tol);
    _, state = iterate(iterator);
    return (iterator, state)
end

function solve!(iterator::Iterator, state::IteratorState)

    while true
        next = iterate(iterator, state);
        
        if next isa Nothing
            break
        else
            _, state = next;
        end
    end
    return state
end

solve(problem::Problem, algorithm::Algorithm; i_max::Integer=10_000, g_tol::Real=1e-8) = solve!(
    init(problem, algorithm; i_max, g_tol)...
);