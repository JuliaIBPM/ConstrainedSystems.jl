# ConstrainedSystems.jl
_Tools for solving constrained dynamical systems_

| Documentation | Build Status |
|:---:|:---:|
| [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaIBPM.github.io/ConstrainedSystems.jl/stable) [![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaIBPM.github.io/ConstrainedSystems.jl/dev) | [![Build Status](https://github.com/JuliaIBPM/ConstrainedSystems.jl/workflows/CI/badge.svg)](https://github.com/JuliaIBPM/ConstrainedSystems.jl/actions) [![Coverage](https://codecov.io/gh/JuliaIBPM/ConstrainedSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIBPM/ConstrainedSystems.jl) |


This package contains several tools for solving and advancing (large-scale) dynamical systems with constraints. These systems generically have the form 

dy/dt = L y - B<sub>1</sub><sup>T</sup> z + r<sub>1</sub>(y,t)

B<sub>2</sub> y = r<sub>2</sub>(y,t)

y(0) = y<sub>0</sub>

where y is a state vector, L is a linear operator with an associated matrix exponential (integrating factor), and z is a constraint force vector (i.e., Lagrange multipliers).

Some of the key components of this package are

* Tools for solving linear algebra problems with constraints and associated Lagrange multipliers, known generically as *saddle point systems*. The sizes of these systems might be large.

* Time integrators that can incorporate these constraints, such as half-explicit Runge-Kutta (HERK) and integrating factor Runge-Kutta (IFRK), or their combination (IF-HERK). These
extend the tools in the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) package, and utilize the same basic syntax for setting
up a problem and solving it.

* Allowance for variable constraint operators B<sub>1</sub><sup>T</sup> and B<sub>2</sub>,
through the use of a variable parameter argument and an associated parameter update
function.

* The ability to add an auxiliary (unconstrained) system of equations that the
constraint operators B<sub>1</sub><sup>T</sup> and B<sub>2</sub> depend upon.

The package is agnostic to the type of systems, and might arise from, e.g., fluid dynamics or rigid-body mechanics.
