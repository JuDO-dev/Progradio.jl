module Progradio

include("base.jl")

abstract type ProgradioProblem{F<:AbstractFloat, I<:Integer} end
abstract type ProgradioDirection{F<:AbstractFloat, I<:Integer} end
abstract type ProgradioSearch{F<:AbstractFloat, I<:Integer} end

abstract type ProgradioDirectionState{F<:AbstractFloat, I<:Integer} end
abstract type ProgradioSearchState{F<:AbstractFloat, I<:Integer} end
abstract type ProgradioIteratorState{F<:AbstractFloat, I<:Integer, DS<:ProgradioDirectionState{F, I}, SS<:ProgradioSearchState{F, I}} end

#include("iterator.jl")

#include("problems/unconstrained.jl")
#include("problems/boxConstrained.jl")
#include("problems/simplexBoxConstrained.jl")

#include("directions/steepestDescent.jl")
#include("directions/conjugateGradient.jl")
#include("directions/quasiNewton.jl")

#include("searches/Armijo.jl")
#include("searches/ArmijoBox.jl")
#include("optimizers/ArmijoSimplexBox.jl")

#include("optimizers/Wolfe.jl")

#include("optimizers/trustRegion.jl")

#include("solve.jl")
#include("optimality.jl")

#include("zoo.jl")

#export UProblem, BCProblem, SBCProblem,
    #SteepestDescent,
    #FletcherReeves, PolakRibiere, HagerZhang,
    #LBFGS,
    #Armijo, #Wolfe, TrustRegion,
    #Iterator#, solve#, optimality, solve_to_optimality
end