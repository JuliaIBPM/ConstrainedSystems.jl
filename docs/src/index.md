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

$$\ddt{y} = L u - B_{1}^{T} z + r_{1}(y,t), \quad B_{2} y = r_{2}(t), \quad y(0) = y_{0}$$

where $y$ is a state vector, $L$ is a linear operator with an associated matrix exponential (integrating factor), and $z$ is a constraint force vector (i.e., Lagrange multipliers). Systems of this type might arise from, e.g., incompressible fluid dynamics, rigid-body mechanics,
or couplings of such systems.

Some of the key components of this package are

* Tools for solving linear algebra problems with constraints and associated Lagrange multipliers, known generically as *saddle point systems*. The sizes of these systems might be large.

* Time integrators that can incorporate these constraints, such as half-explicit Runge-Kutta (HERK) and integrating factor Runge-Kutta (IFRK), or their combination (IF-HERK). These
extend the tools in the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package, and utilize the same basic syntax for setting
up a problem and solving it.

* Allowance for variable constraint operators $B_1^T$ and $B_2$,
through the use of a variable parameter argument and an associated parameter update
function.

* The ability to add an auxiliary (unconstrained) system of equations and state that the
constraint operators $B_1^T$ and $B_2$ and the right-hand side $r_2$ depend upon.


## Installation

This package works on Julia `1.6` and above and is registered in the general Julia registry. To install from the REPL, type
e.g.,
```julia
] add ConstrainedSystems
```

Then, in any version, type
```julia
julia> using ConstrainedSystems
```

The plots in this documentation are generated using [Plots.jl](http://docs.juliaplots.org/latest/). You might want to install that, too, to follow the examples.
