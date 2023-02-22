struct UProblem{F} <: ProgradioProblem{F}
    x_0::Vector{F}  #guess
    f::Function     #objective
    g!::Function    #gradient g!(gx, x)

    function UProblem(x_0::Vector{F}, f::Function, g!::Function) where {F<:AbstractFloat}
        return new{F}(x_0, f, g!)
    end
end

# Pretty printing
Base.show(io::IO, problem::ProgradioProblem) = print(io, typeof(problem), " with ", length(problem.x_0), " variables ");