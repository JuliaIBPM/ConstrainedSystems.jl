# ConstrainedSystems.jl
*Tools for solving constrained dynamical systems*

This package contains several tools for solving and advancing (large-scale) dynamical systems with constraints. Some of the key components of this package are

* Tools for solving linear algebra problems with constraints and associated Lagrange multipliers, known generically as *saddle point systems*. The sizes of these systems might be large.

* Time integrators that can incorporate these constraints, such as half-explicit Runge-Kutta (HERK) and integrating factor HERK.

The package is agnostic to the type of systems, and might arise from, e.g., fluid dynamics or rigid-body mechanics.
