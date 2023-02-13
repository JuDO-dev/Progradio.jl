# Inner product
function dot(a::Vector{F}, b::Vector{F}) where {F<:AbstractFloat}
    dot_sum = zero(F);
    for j in eachindex(a, b)
        dot_sum += a[j] * b[j];
    end
    return dot_sum
end

# Inner product over an index set
function dot(a::Vector{F}, b::Vector{F}, set::BitVector) where {F<:AbstractFloat}
    dot_sum = zero(F);
    for j in eachindex(a, b, set)
        if set[j]
            dot_sum += a[j] * b[j];
        end
    end
    return dot_sum
end

# 2-norm
function norm(a::Vector{F}) where {F<:AbstractFloat}
    sum_of_squares = zero(F);
    for j in eachindex(a)
        sum_of_squares += a[j] ^ 2;
    end
    return sqrt(sum_of_squares)
end

# 2-norm over an index BitVector
function norm(a::Vector{F}, set::BitVector) where {F<:AbstractFloat}
    sum_of_squares = zero(F);
    for j in eachindex(a, set)
        if set[j]
            sum_of_squares += a[j] ^ 2;
        end
    end
    return sqrt(sum_of_squares)
end

# 2-norm over an index BitSet
function norm(a::Vector{F}, set::BitSet) where {F<:AbstractFloat}
    sum_of_squares = zero(F);
    for j in set
        sum_of_squares += a[j] ^ 2;
    end
    return sqrt(sum_of_squares)
end

# Box projection
function project!(x::Vector{F}, ℓ::Vector{F}, u::Vector{F}) where {F<:AbstractFloat}
    for j in eachindex(x, ℓ, u)
        if x[j] <= ℓ[j]
            x[j] = ℓ[j];
        elseif x[j] >= u[j]
            x[j] = u[j];
        end
    end
    return nothing
end

# Box projection over an index set
function project!(x::Vector{F}, ℓ::Vector{F}, u::Vector{F}, set::BitSet) where {F<:AbstractFloat}
    for j in set
        if x[j] <= ℓ[j]
            x[j] = ℓ[j];
        elseif x[j] >= u[j]
            x[j] = u[j];
        end
    end
    return nothing
end

# Orthant projection over an index set
function project!(x::Vector{F}, set::BitSet) where {F<:AbstractFloat}
    for j in set
        if x[j] < 0
            x[j] = zero(F);
        end
    end
    return nothing
end