module ConstrainedSystems

using LinearMaps
#using KrylovKit
using RecursiveArrayTools
#using IterativeSolvers
#using UnPack
using Reexport
@reexport using OrdinaryDiffEq

import OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, alg_order, alg_cache,
                    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                    initialize!, perform_step!, @muladd, @unpack, constvalue,
                    full_cache, @..


import OrdinaryDiffEq.DiffEqBase: AbstractDiffEqLinearOperator,
                                  DEFAULT_UPDATE_FUNC, has_exp,
                                  AbstractODEFunction, isinplace, numargs

import RecursiveArrayTools: recursivecopy, recursivecopy!, recursive_mean


using LinearAlgebra
import LinearAlgebra: ldiv!, mul!, *, \, I

import Base: size, eltype, *, /, +, -

export SaddleSystem, SaddleVector, state, constraint, aux_state, linear_map
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
