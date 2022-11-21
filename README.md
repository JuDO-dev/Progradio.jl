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

## 📦Box-Constrained Problems
For minimizing a function subject to lower and upper bounds on the variables.

$$
\begin{aligned}
\min_x \quad      &f(x)\\
\text{s.t.} \quad &\ell \leq x \leq u,
\end{aligned}
$$

where $x, \ell, u \in \mathbb{R}^n$, and $f: \mathbb{R}^n \rightarrow \mathbb{R}$ is smooth. Defined as:

```julia
bcp = BCProblem(x_0, ℓ, u, f, g!);
```
where `x_0` is a feasible initial guess, and `g!` is the in-place function for the gradient $\nabla f$.

## 📐Simplex-Box-Constrained Problems
For minimizing a function with some of the variables constrained to the unit simplex. The remaining varaibles are constrained by lower and upper bounds.

$$
\begin{aligned}
\min_x \quad        &f(x)\\
\text{s.t.} \quad   &\sum_{j \in \mathcal{S}} x_j = 1, \quad x_j \geq 0 &\forall j \in \mathcal{S},\\
                    &\ell_j \leq x_j \leq u_j &\forall j \notin \mathcal{S},
\end{aligned}
$$

where $\mathcal{S}$ is a set of indices of $x$ in the unit simplex. Defined as:

```julia
sbcp = SBCProblem(x_0, S, ℓ, u, f, g!);
```



## Optimizers

The following methods are implemented for `BCProblem`:

- `Armijo()` is a line-search optimizer with an Armijo-like[^Bertsekas] condition;
- [WIP]`Wolfe()` is a line-search optimizer with a Wolfe-like condition;
- [WIP]`TrustRegion()`.

Each used with a descent direction method:

| Direction | Summary | Variants |
| --- | --- | --- |
| Steepest Descent | Requires `g!` | `SteepestDescent()`
| Conjugate Gradient[^Schwartz] | Restarts every $R$ iterations <br> Requires `g!` | `FletcherReeves()` <br> `PolakRibiere()` <br>  `HagerZhang()`[^Hager] |
| Quasi-Newton | Stores $M$ updates <br> Requires `g!` | [WIP]`LBFGS(M)`[^Nocedal] |

### Recommended usage with `solve(bcp, optimizer)`
```julia
# Problem
bcp = BCProblem(f, g!, x_ℓ, x_u, x_0);

# Solve
solve(bcp, Armijo(FletcherReeves(10)))
```

### Advanced usage with `iterator(bcp, optimizer)`
```julia
# Problem
bcp = BCProblem(f, g!, x_ℓ, x_u, x_0);

# Iterator
bci = iterator(bcp, Armijo(FletcherReeves(10)));

# Iterate
collect(bci)
```

[^Bertsekas]: D. P. Bertsekas, "Projected Newton methods for optimization problems with simple constraints", SIAM Journal on Control and Optimization, Vol. 20, pp.221-246, 1982.

[^Schwartz]: A. Schwartz and E. Polak, "A family of projected descent methods for optimization problems with simple bounds", Journal of Optimization Theory and Applications, Vol. 92, No. 1, pp. 1-32, 1997.

[^Hager]: W. w. Hager and H. Zhang, "A new conjugate gradient method with guaranteed descent and an efficient line search", SIAM Journal on Optimization, Vol. 16, pp. 170-192, 2005.

[^Nocedal]: J. Nocedal, "Updating quasi-Newton matrices with limited storage", Mathematics of Computation, Vol. 35, pp. 773-782, 1980.