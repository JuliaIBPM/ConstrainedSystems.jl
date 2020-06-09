# ConstrainedSystems.jl
_Tools for solving constrained dynamical systems_


| Documentation | Build Status |
|:---:|:---:|
|  | [![Build Status](https://travis-ci.com/JuliaIBPM/ConstrainedSystems.svg?branch=master)](https://travis-ci.com/JuliaIBPM/ConstrainedSystems.jl) [![Build status](https://ci.appveyor.com/api/projects/status/6tokpjqb4x8999g0?svg=true)](https://ci.appveyor.com/project/JuliaIBPM/constrainedsystems-jl) [![codecov](https://codecov.io/gh/JuliaIBPM/ConstrainedSystems.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaIBPM/ConstrainedSystems.jl) |


This package contains several tools for solving and advancing (large-scale) dynamical systems with constraints. Some of the key components of this package are

* Tools for solving linear algebra problems with constraints and associated Lagrange multipliers, known generically as *saddle point systems*. The sizes of these systems might be large.

* Time integrators that can incorporate these constraints, such as half-explicit Runge-Kutta (HERK) and integrating factor HERK.

The package is agnostic to the type of systems, and might arise from, e.g., fluid dynamics or rigid-body mechanics.
