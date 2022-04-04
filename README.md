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

## Problems with Simple Bounds
Also known as box-constrained optimisation, of the form
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\underset{x}{\text{min}} \, f(x) \quad \text{s.t.} \, \, x_{\ell} \leq x \leq x_u" height=24>
</p>

### Solvers
- `ConjugateGradient()`[^Schwartz] uses conjugate directions;
- `L-BFGS()`
- `Newton()`

### Line-search
- `Armijo()` back-tracks until an Armijo-like condition[^Bertsekas] is satisfied. Reliable;
- `Wolfe()` 

### Example
```julia
using Progradio

problem = SBProblem(f, ∇f!, x_0, x_ℓ, x_u);
I = iterator(problem, ConjugateGradient(), maxIterations=12);
collect(I)
```

## Problems with Nonlinear Constraints
<p align="center">
<img src="https://render.githubusercontent.com/render/math?math=\underset{x}{\text{min}} \, f(x) \quad \text{s.t.} \, \, c(x) = 0, x_{\ell} \leq x \leq x_u" height=24>
</p>

## Installation
Using the REPL
```julia
] add Progradio
```


[^Schwartz]: Schwartz, A., and Polak, E., Family of Projected Descent Methods for Optimization Problems with Simple Bounds, Journal of Optimization Theory and Applications, Vol. 92, No. 1, pp. 1-31, 1997. 

[^Bertsekas]: Bertsekas, D. P., Projected Newton Methods for Optimization Problems with Simple Constraints, SIAM Journal on Control and Optimization, Vol. 20,
pp. 221-246, 1982.