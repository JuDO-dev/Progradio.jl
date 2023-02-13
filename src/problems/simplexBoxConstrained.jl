struct SBCProblem{F, I} <: ProgradioProblem{F, I}
    x_0::Vector{F}  #initial guess
    S::BitSet       #simplex indices
    S̄::BitSet       #non-simplex indices
    S_tol::F        #simplex binding tolerance
    ℓ::Vector{F}    #lower bounds
    u::Vector{F}    #upper bounds
    B_tol::F        #box binding tolerance
    f::Function     #objective
    g!::Function    #gradient (mutating function)

    function SBCProblem(x_0::Vector{F}, S::BitSet, S_tol::F, ℓ::Vector{F}, u::Vector{F}, B_tol::F, f::Function, g!::Function; integer_type::Type{I}=Int64) where 
        {F<:AbstractFloat, I<:Integer}
        
        S̄ = BitSet(setdiff(1:length(x_0), S));

        # Ensure that Σx_0_j = 1, x_0_j ≥ 0 ∀j ∈ S
        sum_simplex = zero(F);
        for j in S
            if x_0[j] < 0
                    throw(DomainError(x_0[j], "simplex components of x_0 must be non-negative"))
            end
            sum_simplex += x_0[j];
        end
        if !(sum_simplex ≈ 1.0)
            throw(DomainError("simplex components of x_0 must sum up to 1"))
        end

        # Ensure sizes of x_0, ℓ, and u match
        if !(length(x_0) == length(ℓ) == length(u))
            throw(DimensionMismatch("sizes of x_0, ℓ, and u must match"))
        end

        #Ensure that ℓ_j ≤ x_0_j ≤ u_j ∀j ∉ S
        for j in S̄
            if !(ℓ[j] ≤ x_0[j] ≤ u[j])
                throw(DomainError(x_0[j], "non-simplex components of x_0 must satisfy box constraints"))
            end
        end

        # Ensure that B_tol and S_tol are non-negative
        if B_tol < 0 || S_tol < 0
            throw(DomainError((B_tol, S_tol), "binding tolerances B_tol and S_tol must be positive"))
        end

        # Proceed with construction
        return new{F, I}(x_0, S, S̄, S_tol, ℓ, u, B_tol, f, g!)
    end
end

function SBCProblem(x_0::Vector{F}, S::BitSet, ℓ::Vector{F}, u::Vector{F}, f::Function, g!::Function; integer_type::Type{I}=Int64) where 
    {F<:AbstractFloat, I<:Integer}
    
    # Binding tolerance 10% of smallest bound range
    B_tol = one(F);
    S̄ = BitSet(setdiff(1:length(x_0), S));
    for j in S̄
        B_tol_j = abs(u[j] - ℓ[j]);
        if B_tol_j < B_tol && !isapprox(B_tol_j, 0)
            B_tol = B_tol_j;
        end
    end
    B_tol *= convert(F, 0.1);
    S_tol = convert(F, 0.1);

    return SBCProblem(x_0, S, S_tol, ℓ, u, B_tol, f, g!; integer_type)
end

# Pretty printing
Base.show(io::IO, problem::SBCProblem) = print(io, typeof(problem), " with ", length(problem.x_0), " variables, ", length(problem.S), " of which in the unit simplex ");