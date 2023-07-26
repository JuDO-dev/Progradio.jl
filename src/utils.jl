dot(a::AbstractArray{T, NA}, b::AbstractArray{T, NB}) where {T<:Number, NA, NB} = sum(a[j] * b[j] for j in eachindex(a, b));

dot(a::AbstractArray{T, NA}, b::AbstractArray{T, NB}, bit_vector::BitVector) where {T<:Number, NA, NB} = sum(
    a[j] * b[j] for j in eachindex(a, b, bit_vector) if bit_vector[j];
    init=zero(T)
);

norm2(a::AbstractArray{T, N}) where {T<:Number, N} = sqrt(sum(a[j]^2 for j in eachindex(a)));

norm2(a::AbstractArray{T, N}, bit_vector::BitVector) where {T<:Number, N} = sqrt(sum(
    a[j]^2 for j in eachindex(a, bit_vector) if bit_vector[j];
    init=zero(T)
));

function binding_from_indices!(W::BitVector, κ_index::Union{Vector{I},SubArray{I, 1, Vector{I}, Tuple{UnitRange{Int}},true}}, n::I) where I<:Integer
    W = falses(n);
    for i in κ_index
        if i <= n
            W[i] = true;
        end
    end
    return nothing
end