using Test
using Progradio
const P = Progradio
Float = typeof(1.0);

include("base.jl")

include("problems/unconstrained.jl")
include("problems/boxConstrained.jl")
#include("problems/simplexBoxConstrained.jl")

include("iterator.jl")

include("directions/steepestDescent.jl")