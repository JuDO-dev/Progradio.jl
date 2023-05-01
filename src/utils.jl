dot(a::AbstractVector{<:Real}, b::AbstractVector{<:Real}) = sum(a[j] * b[j] for j in eachindex(a, b));

dot(a::AbstractVector{R}, b::AbstractVector{R}, bit_vector::BitVector) where {R<:Real} = sum(
    a[j] * b[j] for j in eachindex(a, b, bit_vector) if bit_vector[j];
    init=zero(R)
);

norm2(a::AbstractVector{<:Real}) = sqrt(sum(a[j]^2 for j in eachindex(a)));

norm2(a::AbstractVector{R}, bit_vector::BitVector) where {R<:Real} = sqrt(sum(
    a[j]^2 for j in eachindex(a, bit_vector) if bit_vector[j];
    init=zero(R)
));