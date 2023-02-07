module Progradio

include("base.jl")

abstract type ProgradioProblem{F<:AbstractFloat, I<:Integer} end
include("problems/unconstrained.jl")
include("problems/boxConstrained.jl")
include("problems/simplexBoxConstrained.jl")

abstract type ProgradioDirection{F<:AbstractFloat, I<:Integer} end
abstract type ProgradioDirectionState{F<:AbstractFloat, I<:Integer} end

abstract type ProgradioSearch{F<:AbstractFloat, I<:Integer} end
abstract type ProgradioSearchState{F<:AbstractFloat, I<:Integer} end

include("iterator.jl")

include("directions/steepestDescent.jl")
include("directions/conjugateGradient.jl")
#include("directions/quasiNewton.jl")

include("searches/Armijo.jl")
include("searches/ArmijoBox.jl")
include("searches/ArmijoSimplexBox.jl")

#include("optimizers/Wolfe.jl")

#include("optimizers/trustRegion.jl")

#include("solve.jl")
#include("optimality.jl")

#include("zoo.jl")

export UProblem, BCProblem, SBCProblem,
    Iterator,
    SteepestDescent,
    ConjugateGradient, FletcherReeves, PolakRibiere, HagerZhang,
    #QuasiNewton, 
    Armijo#, #Wolfe, TrustRegion,
    #solve, optimality, solve_to_optimality
end