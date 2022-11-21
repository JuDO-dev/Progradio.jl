module Progradio

#using LinearAlgebra: dot, norm

#include("operators.jl")

abstract type ProgradioProblem{F<:AbstractFloat, I<:Integer} end
include("problems/bcp.jl")
include("problems/sbcp.jl")

#abstract type ProgradioDirection end
#abstract type ProgradioDirectionState{F<:AbstractFloat} end
#include("directions/steepestDescent.jl")
#include("directions/conjugateGradient.jl")
#include("directions/quasiNewton.jl")

#abstract type ProgradioOptimizer{F<:AbstractFloat, D<:ProgradioDirection} end
#abstract type ProgradioOptimizerState{F<:AbstractFloat, DS<:ProgradioDirectionState} end
#include("optimizers/Armijo.jl")
#include("optimizers/Wolfe.jl")
#=
include("iterator.jl")
include("solve.jl")
include("optimality.jl")
include("zoo.jl")
=#
export BCProblem, SBCProblem
#=, iterator,
    SteepestDescent,
    FletcherReeves, PolakRibiere, HagerZhang,
    #LBFGS,
    Armijo, #Wolfe, TrustRegion,
    solve#, optimality, solve_to_optimality
=#
end