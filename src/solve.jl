function solve(bcp::BCProblem, optimizer::ProgradioOptimizer; i_max=20, f_tol=1e-6, g_tol=1e-6, x_tol=1e-6)

    # Iterator
    bci = BCIterator(bcp, optimizer, i_max, f_tol, g_tol, x_tol);

    # 0th iteration
    _, state = iterate(bci);

    # Subsequent iterations
    while true
        iterate!(bci, state);
        println(state.fx)

        if state.i >= bci.i_max
            break
        elseif state.i > 1 && converged!(state, f_tol, g_tol, x_tol)
            break
        end
    end

    return state
end