struct BCProblem{F, I} <: ProgradioProblem{F, I}
    x_0::Vector{F}  #initial guess
    ℓ::Vector{F}    #lower bounds
    u::Vector{F}    #upper bounds
    B_tol::F        #binding set tolerance
    f::Function     #objective
    g!::Function    #gradient (mutating function)
    
    function BCProblem(x_0::Vector{F}, ℓ::Vector{F}, u::Vector{F}, B_tol::F, f::Function, g!::Function; integer_type::Type{I}=Int64) where {F<:AbstractFloat, I<:Integer}
        
        # Ensure sizes of x_0, ℓ, and u match
        if !(length(x_0) == length(ℓ) == length(u))
            throw(DimensionMismatch("sizes of x_0, ℓ, and u must match."))
        end
        
        # Ensure that ℓ ≤ x_0 ≤ u
        for j in eachindex(x_0, ℓ, u)
            if !(ℓ[j] ≤ x_0[j] ≤ u[j])
                throw(DomainError(x_0[j], string("x_0[", j, "] must satisfy box constraints.")))
            end
        end

        # Ensure that B_tol is non-negative
        if B_tol < zero(F)
            throw(DomainError(B_tol, "binding tolerance B_tol must be positive"))
        end

        # Proceed with construction
        return new{F, I}(x_0, ℓ, u, B_tol, f, g!)
    end
end

function BCProblem(x_0::Vector{F}, ℓ::Vector{F}, u::Vector{F}, f::Function, g!::Function; integer_type::Type{I}=Int64) where {F<:AbstractFloat, I<:Integer}
    # Binding tolerance 10% of smallest bound range
    B_tol = one(F);
    for j in eachindex(ℓ, u)
        B_tol_j = abs(u[j] - ℓ[j]);
        if B_tol_j < B_tol && !isapprox(B_tol_j, zero(F))
            B_tol = B_tol_j;
        end
    end
    B_tol *= convert(F, 0.1);

    return BCProblem(x_0, ℓ, u, B_tol, f, g!; integer_type)
end