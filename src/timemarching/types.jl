using OrdinaryDiffEq
import OrdinaryDiffEq: OrdinaryDiffEqAlgorithm, alg_order, alg_cache,
                    OrdinaryDiffEqMutableCache, OrdinaryDiffEqConstantCache,
                    initialize!, perform_step!, @muladd, @unpack, constvalue,
                    full_cache, @..

import OrdinaryDiffEq.DiffEqBase: AbstractDiffEqLinearOperator,
                                  DEFAULT_UPDATE_FUNC, has_exp,
                                  AbstractODEFunction, isinplace

export DiffEqLinearOperator, ConstrainedODEFunction

#### Operator and function types ####

mutable struct DiffEqLinearOperator{T,aType} <: AbstractDiffEqLinearOperator{T}
    A :: aType
    DiffEqLinearOperator(A::aType; update_func=DEFAULT_UPDATE_FUNC,
                          dtype=Float64) where {aType} = new{dtype,aType}(A)
end

(f::DiffEqLinearOperator)(du,u,p,t) = mul!(du.x[1],f.A,u.x[1])

import Base: exp
exp(f::DiffEqLinearOperator,args...) = exp(f.A,args...)
has_exp(::DiffEqLinearOperator) = true


"""
        ConstrainedODEFunction(r1,r2,B1,B2[,A=I])

This specifies the functions and operators that comprise an ODE problem with the form

``
\\dfrac{dy}{dt} = Ay - B_1 z + r_1(y,t)
``

``
B_2 y = r_2(t)
``

The optional linear operator `A` defaults to the identity. The `B1` and `B2` functions must be of the respective
in-place forms `B1(dy,z,p)` (to compute the action of `B1` on `z`) and `B2(dz,y,p)` (to compute that action
of `B2` on `y`). The function `r1` must of the in-place form `r1(dy,y,p,t)`, and `r2` must be in the in-place form
`r2(dz,p,t)`.

An optional argument `param_update_func` can be used to set a function that updates problem parameters with
the current solution. This function must take the in-place form `f(q,u,p,t)` to update some `q`. (`q` might
simply be `p`.)
"""
struct ConstrainedODEFunction{iip,static,F1,F2,TMM,C,Ta,Tt,TJ,JVP,VJP,JP,SP,TW,TWt,TPJ,S,TCV,PF} <: AbstractODEFunction{iip}
    odef :: F1
    conf :: F2
    mass_matrix::TMM
    cache::C
    analytic::Ta
    tgrad::Tt
    jac::TJ
    jvp::JVP
    vjp::VJP
    jac_prototype::JP
    sparsity::SP
    Wfact::TW
    Wfact_t::TWt
    paramjac::TPJ
    syms::S
    colorvec::TCV
    param_update_func :: PF
end

function ConstrainedODEFunction(r1,r2,B1,B2,A=I; param_update_func = DEFAULT_UPDATE_FUNC,
                                                 mass_matrix=I,_func_cache=nothing,
                                                 analytic=nothing,
                                                 tgrad = nothing,
                                                 jac = nothing,
                                                 jvp=nothing,
                                                 vjp=nothing,
                                                 jac_prototype = nothing,
                                                 sparsity=jac_prototype,
                                                 Wfact = nothing,
                                                 Wfact_t = nothing,
                                                 paramjac = nothing,
                                                 syms = nothing,
                                                 colorvec = nothing)

    if isinplace(r1,4) && isinplace(r2,3) && isinplace(B1,3) && isinplace(B2,3)
        @inline r1ext(du,u,p,t) = r1(state(du),state(u),p,t)
        @inline r2ext(du,u,p,t) = r2(constraint(du),p,t)
        @inline B1ext_rhs(du,u,p,t) = (dy = state(du); B1(dy,constraint(u),p); dy .*= -1)
        @inline B2ext_rhs(du,u,p,t) = (dz = constraint(du); B2(dz,state(u),p); dz .*= -1)
    elseif isinplace(r1,3) && isinplace(r2,2) && isinplace(B1,2) && isinplace(B2,2)
        @inline r1ext(u,p,t) = r1(state(u),p,t)
        @inline r2ext(u,p,t) = r2(p,t)
        @inline B1ext_rhs(u,p,t) = -B1(constraint(u),p)
        @inline B2ext_rhs(u,p,t) = -B2(state(u),p)
    else
        error("Inconsistent function signatures")
    end

    fill!(_func_cache,0.0)
    odef_nl = SplitFunction(r1ext,B1ext_rhs;_func_cache=deepcopy(_func_cache))
    odef = SplitFunction(A, odef_nl ;_func_cache=deepcopy(_func_cache))
    conf = SplitFunction(r2ext,B2ext_rhs;_func_cache=deepcopy(_func_cache))

    static = param_update_func == DEFAULT_UPDATE_FUNC ? true : false

    ConstrainedODEFunction{isinplace(odef),static,typeof(odef),
                           typeof(conf),typeof(mass_matrix),typeof(_func_cache),
                           typeof(analytic),typeof(tgrad),typeof(jac),typeof(jvp),typeof(vjp),
                           typeof(jac_prototype),typeof(sparsity),
                           typeof(Wfact),typeof(Wfact_t),typeof(paramjac),typeof(syms),
                           typeof(colorvec),typeof(param_update_func)}(odef,conf,mass_matrix,_func_cache,
                           analytic,tgrad,jac,jvp,vjp,jac_prototype,
                           sparsity,Wfact,Wfact_t,paramjac,syms,colorvec,param_update_func)
end

@inline _ode_A!(du,f::ConstrainedODEFunction,u,p,t) = f.odef.f1(du,u,p,t) # A
@inline _ode_r1!(du,f::ConstrainedODEFunction,u,p,t) = f.odef.f2.f.f1(du,u,p,t) # r_1

@inline _ode_neg_B1!(du,f::ConstrainedODEFunction,u,p,t) = f.odef.f2.f.f2(du,u,p,t) # -B_1^T
@inline _constraint_neg_B2!(du,f::ConstrainedODEFunction,u,p,t) = f.conf.f2(du,u,p,t) # -B_2
@inline _constraint_r2!(du,f::ConstrainedODEFunction,u,p,t) = f.conf.f1(du,u,p,t) # r_2


function (f::ConstrainedODEFunction)(du,u,p,t)
    fill!(f.cache,0.0)
    f.odef(f.cache,u,p,t)
    f.conf(du,u,p,t)
    du .+= f.cache
end
