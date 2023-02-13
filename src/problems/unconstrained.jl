struct UProblem{F, I} <: ProgradioProblem{F, I}
    x_0::Vector{F}  #initial guess
    f::Function     #objective
    g!::Function    #gradient (mutating function)

    function UProblem(x_0::Vector{F}, f::Function, g!::Function; integer_type::Type{I}=Int64) where {F<:AbstractFloat, I<:Integer}
        return new{F, I}(x_0, f, g!)
    end
end

# Pretty printing
Base.show(io::IO, problem::ProgradioProblem) = print(io, typeof(problem), " with ", length(problem.x_0), " variables ");