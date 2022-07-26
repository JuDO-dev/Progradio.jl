# 2D Quadratic Problem
function qp2(e::F, x_ℓ::Vector{F}, x_u::Vector{F}, x_0::Vector{F}) where F
    q(x::Vector{F}) = x[1]^2 + e * x[2]^2;
    function ∇q!(∇qx::Vector{F}, x::Vector{F})
        ∇qx[1] = 2 * x[1];
        ∇qx[2] = 2 * e * x[2];
        return nothing
    end
    return BCProblem{F}(q, ∇q!, x_ℓ, x_u, x_0, 0.2)
end

function Rosenbrock2(x_0::Vector{F}, x_ℓ::Vector{F}, x_u::Vector{F}, b::F) where F
    r(x::Vector{F}) = (one(F) - x[1])^2 + b * (x[2] - x[1]^2)^2;
    function ∇r!(∇rx::Vector{F}, x::Vector{F})
        ∇rx[1] = -2 * ((one(F) - x[1]) + 2 * b * x[1] * (x[2] - x[1]^2));
        ∇rx[2] = 2 * b * (x[2] - x[1]^2);
        return nothing
    end
    return BCProblem{F}(r, ∇r!, x_ℓ, x_u, x_0, 0.2)
end