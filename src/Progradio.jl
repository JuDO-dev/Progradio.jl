module Progradio

abstract type Constraints{T<:Real} end
abstract type ConstraintsState end

abstract type Problem{T<:Real, X<:AbstractArray, C<:Constraints} end

abstract type Direction end
abstract type DirectionState{T<:Real} end

abstract type LineSearch end
abstract type LineSearchState{T<:Real} end

abstract type Algorithm end
abstract type AlgorithmState end

include("iterator.jl")
export Iterator

include("constraints/box.jl")
export Box

include("problems/nonlinear.jl")
export NLProblem

include("directions/steepestDescent.jl")
export SteepestDescent
include("directions/conjugateGradient.jl")
export ConjugateGradient, FletcherReeves, PolakRibiere, HagerZhang

include("line_searches/armijo.jl")
export Armijo
include("line_searches/wolfe.jl")
export Wolfe

include("algorithms/two_metric.jl")
export TwoMetric

include("solve.jl")
export solve, init, solve!

include("utils.jl")

end