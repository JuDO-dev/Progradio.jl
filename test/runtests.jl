using Test
using Progradio
const P = Progradio

include("base.jl")

include("problems/unconstrained.jl")
include("problems/boxConstrained.jl")
#include("problems/simplexBoxConstrained.jl")

include("iterator.jl")

include("directions/steepestDescent.jl")