# Component-wise Projection
function projection(x_j::F, x_ℓ_j::F, x_u_j::F) where F<:AbstractFloat
    if x_j <= x_ℓ_j
        return x_ℓ_j
    elseif x_j >= x_u_j
        return x_u_j
    else
        return x_j
    end
end

# Component-wise Binding Set Identification
function binding(x_j::F, x_ℓ_j::F, x_u_j::F, ϵ_j::F, gx_j::F) where F<:AbstractFloat
    if x_ℓ_j <= x_j <= (x_ℓ_j + ϵ_j) && gx_j > 0
        return true
    elseif (x_u_j - ϵ_j) <= x_j <= x_u_j && gx_j < 0
        return true
    else
        return false
    end
end