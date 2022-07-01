abstract type ProgradioProblem{F<:AbstractFloat} end

# Box-Constrained Problem
struct BCProblem{F} <: ProgradioProblem{F}
    f::Function
    g!::Function
    x_ℓ::Vector{F}
    x_u::Vector{F}
    x_0::Vector{F}
    ϵ::F
end

Base.show(io::IO, p::ProgradioProblem) = print(io, typeof(p));
Base.show(io::IO, ::MIME"text/plain", p::ProgradioProblem) =
    print(io, typeof(p), " with ", length(p.x_0), " variables");