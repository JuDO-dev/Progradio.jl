# Component-wise Projection
function Projection(x_j::T, x_ℓ_j::T, x_u_j::T) where T<:AbstractFloat
    if x_j <= x_ℓ_j
        return x_ℓ_j
    elseif x_j >= x_u_j
        return x_u_j
    else
        return x_j
    end
end

# Component-wise Binding Set Identification
function Binding(x_j::T, x_ℓ_j::T, x_u_j::T, ϵ_j::T, ∇fx_j::T) where T<:AbstractFloat
    if x_ℓ_j <= x_j <= (x_ℓ_j + ϵ_j) && ∇fx_j > zero(T) 
        return true
    elseif (x_u_j - ϵ_j) <= x_j <= x_u_j && ∇fx_j < zero(T)
        return true
    else
        return false
    end
end