module Progradio

using LinearAlgebra: dot, norm

abstract type ProgradioDirection{F<:AbstractFloat} end
abstract type ProgradioDirectionState{F<:AbstractFloat} end

abstract type ProgradioOptimizer{F<:AbstractFloat, D<:ProgradioDirection} end
abstract type ProgradioOptimizerState{F<:AbstractFloat, DS<:ProgradioDirectionState} end

include("operators.jl")
include("problems.jl")
include("iterators.jl")

#include("directions/steepestDescent.jl")
#include("directions/conjugateGradient.jl")
#include("directions/quasiNewton.jl")

#include("optimizers/lineSearch.jl")
#include("optimizers/Armijo.jl")
#include("optimizers/Wolfe.jl")
#include("optimizers/trustRegion.jl")

#include("optimality.jl")
#include("solve.jl")

#include("zoo.jl")

export BCProblem, iterator#,
    #SteepestDescent,
    #FletcherReeves, PolakRibiere, HagerZhang,
    #LBFGS,
    #Armijo,
    #optimality, solve, solve_to_optimality
end