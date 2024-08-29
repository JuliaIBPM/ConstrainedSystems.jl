module ConstrainedSystems

using LinearMaps
#using KrylovKit
using RecursiveArrayTools
#using IterativeSolvers
#using UnPack
using Reexport
#@reexport using OrdinaryDiffEq
@reexport using OrdinaryDiffEqCore
@reexport using OrdinaryDiffEqTsit5

import OrdinaryDiffEqCore: OrdinaryDiffEqAlgorithm, alg_order, alg_cache,
                    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                    initialize!, perform_step!, constvalue,
                    full_cache, get_fsalfirstlast, @muladd, @unpack, @..  


import DiffEqBase: AbstractDiffEqLinearOperator,
                                  DEFAULT_UPDATE_FUNC, has_exp,
                                  AbstractODEFunction, isinplace, numargs

import LinearMaps: LinearMap, FunctionMap

import RecursiveArrayTools: recursivecopy, recursivecopy!, recursive_mean


using LinearAlgebra
import LinearAlgebra: ldiv!, mul!, *, \, I

import Base: size, eltype, *, /, +, -

export SaddleSystem, SaddleVector, state, constraint, aux_state, linear_map
export constraint_from_state!
export solvector, mainvector
export SchurSolverType, Direct, CG, GMRES, Iterative


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
