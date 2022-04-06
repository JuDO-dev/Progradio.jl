# Component-wise Projection Operator
function Projection(x_i::T, ℓ_i::T, u_i::T) where {T<:AbstractFloat}
    if x_i <= ℓ_i
        return ℓ_i
    elseif x_i >= u_i
        return u_i
    else
        return x_i
    end
end

# Component-wise Binding Set
function Binding(x_i::T, ℓ_i::T, u_i::T, ϵ_i::T, ∇fx_i::T) where {T<:AbstractFloat}
    if ℓ_i <= x_i <= (ℓ_i + ϵ_i) && ∇fx_i > zero(T) 
        return true
    elseif (u_i - ϵ_i) <= x_i <= u_i && ∇fx_i < zero(T)
        return true
    else
        return false
    end
end