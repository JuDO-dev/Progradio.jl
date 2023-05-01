abstract type ConjugateGradientVariant end

struct ConjugateGradient{F<:AbstractFloat, V<:ConjugateGradientVariant} <: Direction
    restart::Int
    σ_min::F
    σ_max::F
    variant::V

    function ConjugateGradient(restart::Integer, σ_min::F, σ_max::F, variant::V) where {F<:AbstractFloat, V<:ConjugateGradientVariant}
        
        restart > 0     ? nothing : throw(DomainError(restart,  "Ensure"));
        0 < σ_min < 1   ? nothing : throw(DomainError(σ_min,    "Ensure"));
        σ_max > 1       ? nothing : throw(DomainError(σ_max,    "Ensure"));

        return new{F, V}(restart, σ_min, σ_max, variant)
    end
end

abstract type ConjugateDirectionState <: DirectionState end

function direction!(state::IteratorState, cg::ConjugateGradient)
    
    cg_state = state.algorithm_state.direction_state;
    #@. cg_state.d_previous = cg_state.d;

    if cg_state.r ≥ cg.restart || state.i < 1
        @. cg_state.d = -state.gx;
        cg_state.r = 0;

    else
        μ = coefficient!(state, cg.variant);
        μ = isnan(μ) || isinf(μ) ? 0.0 : μ;

        for j in eachindex(cg_state.d, state.set_state.W)
            if state.set_state.W[j]
                cg_state.d[j] = μ * cg_state.d[j] - state.gx[j];
            else
                cg_state.d[j] = -state.gx[j];
            end
        end
            
        gx_W_size = dot(state.gx, state.gx, state.set_state.W);
        too_small = -dot(cg_state.d, state.gx, state.set_state.W) ≤ cg.σ_min * gx_W_size;
        too_large = dot(cg_state.d, cg_state.d, state.set_state.W) ≥ cg.σ_max * gx_W_size;

        if too_small || too_large    
            for j in eachindex(cg_state.d, state.set_state.W)
                if state.set_state.W[j]
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

mutable struct FletcherReevesState{X<:AbstractVector} <: ConjugateDirectionState
    r::Int
    d::X
    gx_previous::X
end

build_state(x::AbstractVector, ::ConjugateGradient{F, V}) where {F, V<:FletcherReeves} = FletcherReevesState(0, zero(x), zero(x));

function coefficient!(state::IteratorState, ::FletcherReeves)

    fr_state = state.algorithm_state.direction_state;
    
    return dot(state.gx, state.gx, state.set_state.W) / dot(fr_state.gx_previous, fr_state.gx_previous, state.set_state.W);
end

struct PolakRibiere <: ConjugateGradientVariant end

mutable struct PolakRibiereState{X<:AbstractVector} <: ConjugateDirectionState
    r::Int
    d::X
    gx_previous::X
    Δx::X
end

build_state(x::AbstractVector, ::ConjugateGradient{F, V}) where {F, V<:PolakRibiere} = PolakRibiereState(0, zero(x), zero(x), zero(x));

function coefficient!(state::IteratorState, ::PolakRibiere)
    
    pr_state = state.algorithm_state.direction_state;
    @. pr_state.Δx = state.gx - pr_state.gx_previous;

    return dot(state.gx, pr_state.Δx, state.set_state.W) / dot(pr_state.gx_previous, pr_state.gx_previous, state.set_state.W);
end

struct HagerZhang <: ConjugateGradientVariant end

mutable struct HagerZhangState{X<:AbstractVector} <: ConjugateDirectionState
    r::Int
    d::X
    d_previous::X
    gx_previous::X
    Δx::X
end

build_state(x::AbstractVector, ::ConjugateGradient{F, V}) where {F, V<:HagerZhang} = HagerZhangState(0, zero(x), zero(x), zero(x), zero(x));

function coefficient!(state::IteratorState, ::HagerZhang)

    hz_state = state.algorithm_state.direction_state;
    @. hz_state.Δx = state.gx - hz_state.gx_previous;
    denominator = dot(hz_state.d_previous, hz_state.Δx, state.set_state.W);
    fraction = dot(hz_state.Δx, hz_state.Δx, state.set_state.W) / denominator;
    @. hz_state.Δx -= 2 * hz_state.d_previous * fraction;

    @. hz_state.d_previous = hz_state.d;
    
    return dot(hz_state.Δx, state.gx, state.set_state.W) / denominator 
end