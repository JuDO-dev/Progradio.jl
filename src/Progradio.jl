module Progradio

using LinearAlgebra: dot, norm

include("problems.jl")
include("operators.jl")

abstract type ProgradioIterator{F<:AbstractFloat} end
Base.length(I::ProgradioIterator) = 1 + I.iterations;
abstract type ProgradioIteratorState{F<:AbstractFloat} end

abstract type ProgradioDirection{F<:AbstractFloat} end
#include("directions/steepestDescent.jl")
#include("directions/conjugateGradient")
#include("directions/quasiNewton.jl")

abstract type ProgradioOptimizer{F<:AbstractFloat} end
#include("optimizers/lineSearch.jl")
#include("optimizers/trustRegion.jl")

#include("optimality.jl")
#include("zoo.jl")

export BCProblem
    #SteepestDescent,    
    #ConjugateGradient, HagerZhang, PolakRibiere, FletcherReeves,
    #QuasiNewton, LBFGS,
    #Armijo,
    #iterator
end