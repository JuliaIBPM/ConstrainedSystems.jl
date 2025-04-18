module ConstrainedSystems

using LinearMaps
using RecursiveArrayTools
using Reexport
using MuladdMacro
@reexport using SciMLBase
using OrdinaryDiffEqCore


import MuladdMacro: @muladd     
import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache, 
        initialize!, perform_step!, @unpack, constvalue, full_cache, @..,
        alg_order,alg_cache,full_cache,get_fsalfirstlast               

using SciMLBase: AbstractSciMLOperator,ODEFunction, AbstractODEFunction, SplitFunction,
                 ODEProblem, init, solve

using DiffEqBase: DEFAULT_UPDATE_FUNC, isinplace, numargs
import DiffEqBase: has_exp
import DiffEqBase: UNITLESS_ABS2, ODE_DEFAULT_NORM, recursive_length

import LinearMaps: LinearMap, FunctionMap

import RecursiveArrayTools: recursivecopy, recursivecopy!, recursive_mean


using LinearAlgebra
import LinearAlgebra: ldiv!, mul!, *, \, I

import Base: size, eltype, *, /, +, -

export SaddleSystem, SaddleVector, state, constraint, aux_state, linear_map
export constraint_from_state!
export solvector, mainvector
export SchurSolverType, Direct, CG, GMRES, Iterative
export DiffEqLinearOperator, ConstrainedODEFunction


include("vectors.jl")
include("saddlepoint/saddlesystems.jl")
include("saddlepoint/linearmaps.jl")
include("saddlepoint/arithmetic.jl")
include("saddlepoint/testproblems.jl")


include("timemarching/types.jl")
include("timemarching/misc_utils.jl")
include("timemarching/timesaddlesystems.jl")
include("timemarching/algorithms.jl")
include("timemarching/testproblems.jl")


end
