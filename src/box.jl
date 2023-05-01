struct Box{R, X} <: ConstraintSet{R, X}
    ℓ::X
    u::X
    B_tol::R

    function Box(ℓ::X, u::X, B_tol::R) where {R<:Real, X<:AbstractVector{R}}

        B_tol ≥ 0               ? nothing : throw(DomainError(B_tol,        "Ensure"));
        length(ℓ) == length(u)  ? nothing : throw(DimensionMismatch(        "Ensure"));
        
        for j in eachindex(ℓ, u)
            u[j] ≥ ℓ[j]         ? nothing : throw(DomainError((ℓ[j], u[j]), "Ensure "));
        end

        return new{R, X}(ℓ, u, B_tol)
    end
end

function Box(ℓ::X, u::X) where {R<:Real, X<:AbstractVector{R}}

    Δx_min = one(R);
    Δx_j = zero(R);
    
    for j in eachindex(ℓ, u)
        Δx_j = abs(u[j] - ℓ[j]);
        Δx_j ≈ 0 ? continue : nothing;
        
        if Δx_j < Δx_min
            Δx_min = Δx_j;
        end
    end

    return Box(ℓ, u, Δx_min / 10)
end

function is_inside(x::X, box::Box{R, X}) where {R<:Real, X<:AbstractVector{R}}

    for j in eachindex(box.ℓ, x, box.u)
        if !(box.ℓ[j] ≤ x[j] ≤ box.u[j])
            return false
        end
    end
    
    return true
end

mutable struct BoxState <: ConstraintSetState
    B_ℓ::BitVector
    B_u::BitVector
    B::BitVector
    W::BitVector
end

build_state(x::AbstractVector, ::Box) = BoxState(falses(length(x)), falses(length(x)), falses(length(x)), trues(length(x)));

function binding!(state::IteratorState{R, X, AS, SS}, box::Box{R, X}) where {R, X, AS, SS}

    box_state = state.set_state;
    
    for j in eachindex(box.ℓ, state.x, box.u)
        if box.ℓ[j] ≤ state.x[j] ≤ (box.ℓ[j] + box.B_tol) && state.gx[j] > 0
            box_state.B_ℓ[j] = true;
            box_state.B_u[j] = false;
            box_state.B[j] = true;
            box_state.W[j] = false;
        
        elseif (box.u[j] - box.B_tol) ≤ state.x[j] ≤ box.u[j] && state.gx[j] < 0
            box_state.B_ℓ[j] = false;
            box_state.B_u[j] = true;
            box_state.B[j] = true;
            box_state.W[j] = false;

        else
            box_state.B_ℓ[j] = false;
            box_state.B_u[j] = false;
            box_state.B[j] = false;
            box_state.W[j] = true;
        end
    end
    
    return nothing
end

function has_converged(state::IteratorState, g_tol::Real, ::Box)

    for j in eachindex(state.gx)
        if state.set_state.W[j]
            if abs(state.gx[j]) > g_tol
                return false
            end
        end
    end
    return true
end

function project!(x::X, box::Box{R, X}) where {R, X}
    
    for j in eachindex(x, box.ℓ, box.u)
        if x[j] <= box.ℓ[j]
            x[j] = box.ℓ[j];
        elseif x[j] >= box.u[j]
            x[j] = box.u[j];
        end
    end

    return nothing
end

function project!(x::X, box::Box{R, X}, bit_set::BitSet) where {R, X}

    for j in bit_set
        if x[j] <= box.ℓ[j]
            x[j] = box.ℓ[j];
        elseif x[j] >= box.u[j]
            x[j] = box.u[j];
        end
    end

    return nothing
end