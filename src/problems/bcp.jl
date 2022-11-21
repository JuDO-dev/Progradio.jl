# Box-constrained problem
struct BCProblem{F, I} <: ProgradioProblem{F, I}
    x_0::Vector{F}  #initial guess
    x_ℓ::Vector{F}  #lower bounds
    x_u::Vector{F}  #upper bounds
    f::Function     #objective
    g!::Function    #gradient (mutating function)
    #H!             #Hessian
    ϵ_B::F          #binding set tolerance
    n_x::I          #number of variables
    
    function BCProblem(x_0::Vector{F}, x_ℓ::Vector{F}, x_u::Vector{F}, f::Function, g!::Function; ϵ_B::F=0.1, integer_type::Type{I}=Int64) where {F<:AbstractFloat, I<:Integer}
        n_x_0 = length(x_0)

        # Check sizes of x_0, x_ℓ, and x_u
        if !(n_x_0 == length(x_ℓ) == length(x_u))
            return error("Sizes of x_0, x_ℓ, and x_u must match")
        
        else
            # Check if x_ℓ ≤ x_0 ≤ x_u
            out_of_box = false;
            for j in 1:n_x_0
                if !(x_ℓ[j] ≤ x_0[j] ≤ x_u[j])
                    out_of_box = true;
                    break
                end
            end
            if out_of_box
                return error("x_ℓ ≤ x_0 ≤ x_u must hold")

            # Proceed with construction
            else
                return new{F, I}(x_0, x_ℓ, x_u, f, g!, ϵ_B, convert(integer_type, n_x_0))
            end
        end
    end
end

# Pretty printing
Base.show(io::IO, bcp::BCProblem) = print(io, typeof(bcp), " with ", bcp.n_x, " variables");