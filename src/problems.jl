struct NLProblem{R, X, S<:ConstraintSet{R, X}} <: Problem{R, X}
    x_guess::X
    set::S
    f::Function
    g!::Function

    function NLProblem(x_guess::X, set::S, f::Function, g!::Function) where {R, X<:AbstractVector{R}, S<:ConstraintSet{R}}

        is_inside(x_guess, set) ? nothing : throw(DomainError(""));

        return new{R, X, S}(x_guess, set, f, g!)
    end
end