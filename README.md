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

## ♾️Unconstrained Problems

$$
\begin{aligned}
\min_x \hspace{0.5em} f(x)
\end{aligned}
$$

where $x \in \mathbb{R}^n$, and $f: \mathbb{R}^n \rightarrow \mathbb{R}$ is smooth. Given an initial guess `x_0::Vector` and an in-place gradient function `g!`, the problem is defined as:
```julia
up = UProblem(x_0, f, g!);
```

## 📦Box-Constrained Problems

$$
\begin{aligned}
\min_x \hspace{0.5em}      &f(x)\\
\text{s.t.} \hspace{0.5em} &\ell \leq x \leq u,
\end{aligned}
$$

where $\ell, u \in \mathbb{R}^n$. Given `ℓ::Vector` and `u::Vector`, the problem is defined as:
```julia
bcp = BCProblem(x_0, ℓ, u, f, g!);
```

## 📐Simplex-Box-Constrained Problems

$$
\begin{aligned}
\min_x \hspace{0.5em}       &f(x)\\
\text{s.t.} \hspace{0.5em}  &\sum_{j \in \mathcal{S}} x_j = 1, \quad x_j \geq 0 &\forall j \in \mathcal{S},\\
                    &\ell_j \leq x_j \leq u_j &\forall j \notin \mathcal{S},
\end{aligned}
$$

where $\mathcal{S}$ is the set of indices of $x$ in the unit simplex. Given `S::BitSet`, the problem is defined as:
```julia
sbcp = SBCProblem(x_0, S, ℓ, u, f, g!);
```

## Available Methods

|Direction\Search|`Armijo()`|`Wolfe()`|`TrustRegion()`|
|:-:|:-:|:-:|:-:|
|`SteepestDescent()`|♾️📦📐[^Bertsekas]|-|-|      
|`ConjugateGradient()`|♾️📦|-|-|
|`QuasiNewton()`|-|-|-|
|`Newton()`|-|-|-|

## Usage
Recommended usage with `solve()`
```julia
# Problem
bcp = BCProblem(x_0, ℓ, u, f, g!);

# Solve
solve(bcp, SteepestDescent(), Armijo())
```
Advanced usage with `Iterator`
```julia
# Iterator
iterator = Iterator(bcp, SteepestDescent(), Armijo());

# Iterate
collect(iterator)
```

[^Bertsekas]: D. P. Bertsekas, "Projected Newton methods for optimization problems with simple constraints", SIAM Journal on Control and Optimization, Vol. 20, pp.221-246, 1982.