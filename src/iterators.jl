abstract type ProgradioIterator{F<:AbstractFloat} end
Base.length(I::ProgradioIterator) = 1 + I.i_max;
Base.eltype(::ProgradioIterator{F}) where F = F;

# Box-constrained problems
struct BCIterator{F, O<:ProgradioOptimizer} <: ProgradioIterator{F}
    # Main
    bcp::BCProblem
    optimizer::O
    # Termination
    i_max::Integer
    f_tol::F
    g_tol::F
    x_tol::F
end

function iterator(bcp::BCProblem{F}, optimizer::O; i_max=20, f_tol=1e-6, g_tol=1e-6, x_tol=1e-6) where {F<:AbstractFloat, O<:ProgradioOptimizer}
    return BCIterator{F, O}(bcp, optimizer, i_max, f_tol, g_tol, x_tol)
end