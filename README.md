[//]: Logo
<p align="center">
    <img src="./docs/src/assets/logo256px.svg">
</p>

# Projected Gradient Optimization

[//]: Badges
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuDO-dev.github.io/Progradio.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuDO-dev.github.io/Progradio.jl/dev)
[![Build Status](https://github.com/JuDO-dev/Progradio.jl/actions/workflows/CI.yml/badge.svg?branch=dev)](https://github.com/JuDO-dev/Progradio.jl/actions/workflows/CI.yml?query=branch%3Adev)
[![Coverage](https://codecov.io/gh/JuDO-dev/Progradio.jl/branch/dev/graph/badge.svg)](https://codecov.io/gh/JuDO-dev/Progradio.jl)

## Installation
```julia
using Pkg; Pkg.add("Progradio")
```

## Box-Constrained Optimization

The problem can be formulated as:
$$\min_x {f(x)} \quad \text{s.t.} \quad x_{\ell} \leq x \leq x_u,$$
where $x, x_{\ell}, x_u \in \mathbb{R}^n$ and $f: \mathbb{R}^n \rightarrow \mathbb{R}$, near an initial guess $x_0$. Implemented as:

```julia
BCProblem(f, g!, x_ℓ, x_u, x_0)
```
where `g!` is the gradient $\nabla f$ defined in-place.  


### Line-Search Optimizer: `Armijo(direction)`[^Bertsekas]

| Direction | Summary | Variants |
| --- | --- | --- |
| Steepest Descent | Requires `g!` | `SteepestDescent()`
| Conjugate Gradient[^Schwartz] | Requires `g!` <br> Restarts every `r` iterations | `HagerZhang(r)`[^Hager] <br> `PolakRibiere(r)` <br> `FletcherReeves(r)` |
| Quasi-Newton | Requires `g!` <br> Stores `m` updates | `LBFGS(m)`[^Nocedal] |

### Recommended usage with `solve(bcp, optimizer)`
```julia
# Problem
bcp = BCProblem(f, g!, x_ℓ, x_u, x_0);

# Solve
solve(bcp, Armijo(HagerZhang(10)))
```

### Advanced usage with `iterator(bcp, optimizer)`
```julia
# Problem
bcp = BCProblem(f, g!, x_ℓ, x_u, x_0);

# Iterator
I = iterator(bcp, Armijo(HagerZhang(10)));

# Iterate
collect(I)
```

[^Bertsekas]: D. P. Bertsekas, "Projected Newton methods for optimization problems with simple constraints", SIAM Journal on Control and Optimization, Vol. 20, pp.221-246, 1982.

[^Schwartz]: A. Schwartz and E. Polak, "A family of projected descent methods for optimization problems with simple bounds", Journal of Optimization Theory and Applications, Vol. 92, No. 1, pp. 1-32, 1997.

[^Hager]: W. w. Hager and H. Zhang, "A new conjugate gradient method with guaranteed descent and an efficient line search", SIAM Journal on Optimization, Vol. 16, pp. 170-192, 2005.

[^Nocedal]: J. Nocedal, "Updating quasi-Newton matrices with limited storage", Mathematics of Computation, Vol. 35, pp. 773-782, 1980.