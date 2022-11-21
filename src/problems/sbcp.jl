# Simplex-box-constrained problem
struct SBCProblem{F, I} <: ProgradioProblem{F, I}
    x_0::Vector{F}  #initial guess
    S::BitVector    #simplex indices
    x_ℓ::Vector{F}  #lower bounds
    x_u::Vector{F}  #upper bounds
    f::Function     #objective
    g!::Function    #gradient (mutating function)
    #H!             #Hessian
    ϵ_B::F          #binding set tolerance
    n_S::I          #number of simplex variables
    n_x::I          #number of variables

    function SBCProblem(x_0::Vector{F}, S::BitVector, x_ℓ::Vector{F}, x_u::Vector{F}, f::Function, g!::Function; ϵ_B::F=0.1, integer_type::Type{I}=Int64) where {F<:AbstractFloat, I<:Integer}
        n_S = sum(S);
        n_x_0 = length(x_0);
        
        # Check number of simplex variables
        if !(2 ≤ n_S ≤ n_x_0)
            return error("Number of simplex components must be between 2 and n_x")

        # Check sizes of x_0, x_ℓ, x_u and S
        elseif !(n_x_0 == length(x_ℓ) == length(x_u) == length(S))
            return error("Sizes of x_0, x_ℓ, x_u, and S must match")
         
        # Check if Σx_0_j = 1, x_0_j ≥ 0 ∀j ∈ S
        # And if x_ℓ_j ≤ x_0_j ≤ x_u_j ∀j ∉ S
        else
            sum_simplex = 0.0;
            out_of_simplex = false;
            out_of_box = false;
            for j in 1:n_x_0
                if S[j]
                    sum_simplex += x_0[j];
                    if x_0[j] < 0.0
                        out_of_simplex = true;
                        break
                    end
                else
                    if !(x_ℓ[j] ≤ x_0[j] ≤ x_u[j])
                        out_of_box = true;
                    end
                end
            end
            if out_of_simplex || !(sum_simplex ≈ 1.0)
                return error("Simplex components of x_0 must be non-negative and sum up to 1.")
            elseif out_of_box
                return error("Non-simplex components of x_0 must satisfy box constraints.")

            # Proceed with construction
            else
                return new{F, I}(x_0, S, x_ℓ, x_u, f, g!, ϵ_B, convert(integer_type, n_S), convert(integer_type, n_x_0))
            end
        end
    end
end

# Pretty printing
Base.show(io::IO, sbcp::SBCProblem) = print(io, typeof(sbcp), " with ", sbcp.n_x, " variables, ", sbcp.n_S, " of which on a simplex");