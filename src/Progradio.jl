module Progradio

# Common
include("operators.jl")
#include("optimality.jl")

# Optimisers
#include("conjugateGradient.jl")
#include("l-BFGS.jl")

# Line-search
#include("Armijo.jl")
#include("Wolfe.jl")

# Interfaces
#include("iterator.jl")
#include("solve.jl")


end