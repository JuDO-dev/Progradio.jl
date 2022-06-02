# Simple Bounds Problem
struct SBProblem{T<:AbstractFloat}
    f
    ∇f!
    x_0::Vector{T}
    x_ℓ::Vector{T}
    x_u::Vector{T}
end

# 2D Quadratic Example
function QP2(x_0::Vector{T}, x_ℓ::Vector{T}, x_u::Vector{T}, c::T) where T<:AbstractFloat
    q(x::Vector{T}) = x[1]^2 + c * x[2]^2;
    function ∇q!(∇qx::Vector{T}, x::Vector{T})
        ∇qx[1] = 2 * x[1];
        ∇qx[2] = 2 * c * x[2];
        return nothing
    end
    return SBProblem{T}(q, ∇q!, x_0, x_ℓ, x_u)
end

# 2D Rosenbrock Example
function Rosenbrock2(x_0::Vector{T}, x_ℓ::Vector{T}, x_u::Vector{T}, b::T) where T<:AbstractFloat
    r(x::Vector{T}) = (one(T) - x[1])^2 + b * (x[2] - x[1]^2)^2;
    function ∇r!(∇rx::Vector{T}, x::Vector{T})
        ∇rx[1] = -2 * ((one(T) - x[1]) + 2 * b * x[1] * (x[2] - x[1]^2));
        ∇rx[2] = 2 * b * (x[2] - x[1]^2);
        return nothing
    end
    return SBProblem{T}(r, ∇r!, x_0, ℓ, u)
end