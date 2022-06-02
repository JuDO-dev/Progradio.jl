module Progradio

using LinearAlgebra: dot, norm

# Common
include("operators.jl")
include("problems.jl")

# Iterator interface
#abstract type ProgradioIterator{T<:AbstractFloat} end
#abstract type ProgradioIteratorState{T<:AbstractFloat} end
#Base.eltype(::ProgradioIterator{T}) where T<:AbstractFloat = T;
#Base.length(I::ProgradioIterator) = 1 + I.maxIterations;

# Line-search
#abstract type ProgradioLineSearch{T<:AbstractFloat} end
#include("Armijo.jl")
#include("Wolfe.jl")

# Optimisers
#abstract type ProgradioOptimiser{T<:AbstractFloat} end
#include("conjugateGradient.jl")
#include("lBFGS.jl")

# Interfaces
#include("solve.jl")

# Utils
#include("optimality.jl")

export SBProblem#,
#    Armijo,
#    ConjugateGradient, FletcherReeves, PolakRibiere

end