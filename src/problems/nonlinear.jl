struct NLProblem{T, X, C} <: Problem{T, X, C}
    x_start::X
    constraints::C
    f::Function
    g!::Function

    NLProblem(x::X, c::C, f, g!) where {T<:Real, N, X<:AbstractArray{T, N},
        C<:Constraints{T}} = new{T, X, C}(x, c, f, g!
    );
end