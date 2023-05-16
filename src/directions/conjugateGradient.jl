abstract type ConjugateGradientVariant end

struct ConjugateGradient{T<:Real, V<:ConjugateGradientVariant} <: Direction
    restart::Int
    σ_min::T
    σ_max::T
    variant::V

    function ConjugateGradient(restart::Integer, σ_min::T, σ_max::T, variant::V) where {T<:Real, V<:ConjugateGradientVariant}
        
        restart > 0     ? nothing : throw(DomainError(restart,  "Ensure"));
        0 < σ_min < 1   ? nothing : throw(DomainError(σ_min,    "Ensure"));
        σ_max > 1       ? nothing : throw(DomainError(σ_max,    "Ensure"));

        return new{T, V}(restart, σ_min, σ_max, variant)
    end
end

abstract type ConjugateDirectionState{T} <: DirectionState{T} end

function direction!(state::IteratorState, cg::ConjugateGradient)
    
    cg_state = state.algorithm_state.direction_state;

    if cg_state.r ≥ cg.restart || state.i < 1
        @. cg_state.d = -state.gx;
        cg_state.r = 0;

    else
        μ = coefficient!(state, cg.variant);
        μ = isnan(μ) || isinf(μ) ? 0.0 : μ;

        for j in eachindex(cg_state.d, state.constraints_state.W)
            if state.constraints_state.W[j]
                cg_state.d[j] = μ * cg_state.d[j] - state.gx[j];
            else
                cg_state.d[j] = -state.gx[j];
            end
        end
            
        gx_W_size = dot(state.gx, state.gx, state.constraints_state.W);
        too_small = -dot(cg_state.d, state.gx, state.constraints_state.W) ≤ cg.σ_min * gx_W_size;
        too_large = dot(cg_state.d, cg_state.d, state.constraints_state.W) ≥ cg.σ_max * gx_W_size;

        if too_small || too_large    
            for j in eachindex(cg_state.d, state.constraints_state.W)
                if state.constraints_state.W[j]
                    cg_state.d[j] = -state.gx[j];
                end
            end
            cg_state.r = 0;
        else
            cg_state.r += 1;
        end
    end

    @. cg_state.gx_previous = state.gx;

    return nothing
end

struct FletcherReeves <: ConjugateGradientVariant end

mutable struct FletcherReevesState{T} <: ConjugateDirectionState{T}
    r::Int
    d::Vector{T}
    gx_previous::Vector{T}
end

build_state(x::X, ::ConjugateGradient{T, V}) where {T<:Real, N, X<:AbstractArray{T, N}, V<:FletcherReeves} = 
    FletcherReevesState(0, Vector{T}(undef, length(x)), Vector{T}(undef, length(x))
);

function coefficient!(state::IteratorState, ::FletcherReeves)

    fr_state = state.algorithm_state.direction_state;
    
    return dot(state.gx, state.gx, state.constraints_state.W) / dot(fr_state.gx_previous, fr_state.gx_previous, state.constraints_state.W);
end

struct PolakRibiere <: ConjugateGradientVariant end

mutable struct PolakRibiereState{T} <: ConjugateDirectionState{T}
    r::Int
    d::Vector{T}
    gx_previous::Vector{T}
    Δgx::Vector{T}
end

build_state(x::X, ::ConjugateGradient{T, V}) where {T<:Real, N, X<:AbstractArray{T, N}, V<:PolakRibiere} =
    PolakRibiereState(0, Vector{T}(undef, length(x)), Vector{T}(undef, length(x)), Vector{T}(undef, length(x))
);

function coefficient!(state::IteratorState, ::PolakRibiere)
    
    pr_state = state.algorithm_state.direction_state;
    @. pr_state.Δgx = state.gx - pr_state.gx_previous;

    return dot(state.gx, pr_state.Δgx, state.constraints_state.W) / dot(pr_state.gx_previous, pr_state.gx_previous, state.constraints_state.W);
end

struct HagerZhang <: ConjugateGradientVariant end

mutable struct HagerZhangState{T} <: ConjugateDirectionState{T}
    r::Int
    d::Vector{T}
    d_previous::Vector{T}
    gx_previous::Vector{T}
    Δgx::Vector{T}
end

build_state(x::X, ::ConjugateGradient{T, V}) where {T<:Real, N, X<:AbstractArray{T, N}, V<:HagerZhang} =
    HagerZhangState(0, Vector{T}(undef, length(x)), Vector{T}(undef, length(x)),
    Vector{T}(undef, length(x)), Vector{T}(undef, length(x))
);

function coefficient!(state::IteratorState, ::HagerZhang)

    hz_state = state.algorithm_state.direction_state;
    @. hz_state.Δgx = state.gx - hz_state.gx_previous;
    denominator = dot(hz_state.d_previous, hz_state.Δgx, state.constraints_state.W);
    fraction = dot(hz_state.Δgx, hz_state.Δgx, state.constraints_state.W) / denominator;
    @. hz_state.Δgx -= 2 * hz_state.d_previous * fraction;

    @. hz_state.d_previous = hz_state.d;
    
    return dot(hz_state.Δgx, state.gx, state.constraints_state.W) / denominator 
end