module Progradio

using LinearAlgebra: dot, norm

abstract type ProgradioDirection{F<:AbstractFloat} end
abstract type ProgradioDirectionState{F<:AbstractFloat} end

abstract type ProgradioOptimizer{F<:AbstractFloat, D<:ProgradioDirection} end
abstract type ProgradioOptimizerState{F<:AbstractFloat, DS<:ProgradioDirectionState} end

include("operators.jl")
include("problems.jl")
include("iterators.jl")

include("directions/steepestDescent.jl")
include("directions/conjugateGradient.jl")
#include("directions/quasiNewton.jl")

include("optimizers/Armijo.jl")
#include("optimizers/Wolfe.jl")
#include("optimizers/trustRegion.jl")

include("solve.jl")
#include("optimality.jl")
include("zoo.jl")

export BCProblem, iterator,
    SteepestDescent,
    FletcherReeves, PolakRibiere, HagerZhang,
    #LBFGS,
    Armijo, #Wolfe, TrustRegion,
    solve#, optimality, solve_to_optimality
end