# ConstrainedSystems.jl

*tools for solving constrained dynamical systems*

```math
\def\ddt#1{\frac{\mathrm{d}#1}{\mathrm{d}t}}

\renewcommand{\vec}{\boldsymbol}
\newcommand{\uvec}[1]{\vec{\hat{#1}}}
\newcommand{\utangent}{\uvec{\tau}}
\newcommand{\unormal}{\uvec{n}}

\renewcommand{\d}{\,\mathrm{d}}
```

This package contains several tools for solving and advancing (large-scale) dynamical systems with constraints. These systems generically have the form

$$\ddt{u} = A u - B_{1}^{T} f + r_{1}(u,t), \quad B_{2} u = r_{2}(u,t), \quad u(0) = u_{0}$$

where $u$ is a state vector, $A$ is a linear operator with an associated matrix exponential (integrating factor), and $f$ is a constraint force vector (i.e., Lagrange multipliers).

Some of the key components of this package are

- Tools for solving linear algebra problems with constraints and associated Lagrange multipliers, known generically as saddle point systems. The sizes of these systems might be large.

- Time integrators that can incorporate these constraints, such as half-explicit Runge-Kutta (HERK) and integrating factor HERK.

The package is agnostic to the type of systems, and might arise from, e.g., fluid dynamics or rigid-body mechanics.

## Installation

This package works on Julia `1.0` and above and is registered in the general Julia registry. To install from the REPL, type
e.g.,
```julia
] add ConstrainedSystems
```

Then, in any version, type
```julia
julia> using ConstrainedSystems
```

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/).
You might want to install that, too, to follow the examples.
