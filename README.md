[//]: Logo
<p align="center">
<img
    src="./docs/src/assets/logo256px.svg"
    width=256px
    >
</p>

# Projected Gradient Optimisation
[//]: Badges
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuDO-dev.github.io/Progradio.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuDO-dev.github.io/Progradio.jl/dev)
[![Build Status](https://github.com/JuDO-dev/Progradio.jl/actions/workflows/CI.yml/badge.svg?branch=dev)](https://github.com/JuDO-dev/Progradio.jl/actions/workflows/CI.yml?query=branch%3Adev)
[![Coverage](https://codecov.io/gh/JuDO-dev/Progradio.jl/branch/dev/graph/badge.svg)](https://codecov.io/gh/JuDO-dev/Progradio.jl)

Solve nonlinear optimisation problems subject to simple bounds (box-constraints), of the form

$$\min_x {f(x)} \quad \text{s.t.} \, \, x_{\ell} \leq x \leq x_u,$$

near an initial guess $x_0$.

## Installation
```julia
using Pkg; Pkg.add("Progradio")
```

## Algorithms

| Type | Summary | Flavours |
| --- | --- | --- |
| `ConjugateGradient()` | Uses the conjugacy of directions; <br> Requires $\nabla f$; <br> Stores $O(n)$ values. | `HagerZhang` <br> `PolakRibiere`[^Schwartz] <br> `FletcherReeves` |
| `QuasiNewton()` | Approximates the Hessian of $f$; <br> Requires $\nabla f$; <br> Stores $O(m \times n)$ values. | `lBFGS`[WIP] |

where $n$ is the problem dimension ($x \in \mathbb{R}^n$) and $m$ is number of stored updates.

### Line-search

- `Armijo()` back-tracks until an Armijo-like condition[^Bertsekas] is satisfied;
- `Wolfe()`[WIP]


## Recommended usage with `solve()`
```julia
using GalacticOptim
using Progradio

```

## Advanced usage with `iterator()`
```julia
using Progradio

# Simple-bounds problem
sbp = SBProblem(f, ∇f!, x_0, x_ℓ, x_u);

# Algorithm
cg = ConjugateGradient(PolakRibiere(), Armijo());

# Iterator
I = iterator(sbp, cg, 20);
collect(I)
```

[^Schwartz]: Schwartz, A., and Polak, E., Family of Projected Descent Methods for Optimization Problems with Simple Bounds, Journal of Optimization Theory and Applications, Vol. 92, No. 1, pp. 1-31, 1997. 

[^Bertsekas]: Bertsekas, D. P., Projected Newton Methods for Optimization Problems with Simple Constraints, SIAM Journal on Control and Optimization, Vol. 20,
pp. 221-246, 1982.