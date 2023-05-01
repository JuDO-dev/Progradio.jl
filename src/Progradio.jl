module Progradio

abstract type ConstraintSet{R<:Real, X<:AbstractVector{R}} end
abstract type ConstraintSetState end

abstract type Problem{R<:Real, X<:AbstractVector{R}} end

abstract type Direction end
abstract type DirectionState end

abstract type LineSearch end
abstract type LineSearchState end

abstract type Algorithm end
abstract type AlgorithmState end

include("iterator.jl")
export Iterator

include("box.jl")
export Box

include("problems.jl")
export NLProblem

include("directions/steepestDescent.jl")
include("directions/conjugateGradient.jl")
export SteepestDescent, ConjugateGradient, FletcherReeves, PolakRibiere, HagerZhang

include("line_searches/armijo.jl")
# # #include("line_searches/Wolfe.jl")
export Armijo#, Wolfe

include("algorithms/two_metric.jl")
export TwoMetric

include("solve.jl")
export solve, init, solve!

include("utils.jl")

end