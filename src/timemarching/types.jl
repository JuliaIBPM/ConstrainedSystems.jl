export DiffEqLinearOperator, ConstrainedODEFunction


#### Operator and function types ####

mutable struct DiffEqLinearOperator{T,aType} <: AbstractDiffEqLinearOperator{T}
    L :: aType
    DiffEqLinearOperator(L::aType; update_func=DEFAULT_UPDATE_FUNC,
                          dtype=Float64) where {aType} = new{dtype,aType}(L)
end

(f::DiffEqLinearOperator)(du,u,p,t) = mul!(du.x[1],f.L,u.x[1])

import Base: exp
exp(f::DiffEqLinearOperator,args...) = exp(f.L,args...)
exp(f::ArrayPartition,t,x::ArrayPartition) =
                ArrayPartition((exp(Li,t,xi) for (Li,xi) in zip(f.L.x,x.x))...)


has_exp(::DiffEqLinearOperator) = true

exp(L::AbstractMatrix,t,x) = exp(factorize(L)*t)
exp(L::UniformScaling,t,x) = exp(Diagonal(L,length(x))*t)

#=
A function of this type should be able to take arguments (du,u,p,t) (for in-place)
or (u,p,t) for out-of-place, and distribute the parts of u and du as
needed to the component functions. u and du will both be of type ArrayPartition,
with two parts: state(u) and constraint(u).

In some cases we wish to solve a separate (unconstrained) system, e.g., to
update the constraint operators B1 and B2. This update is done via the
parameters, using param_update_func. In this case, u and du are ArrayPartition
with the first part corresponding to the constrained system (and of ArrayPartition type)
and the second part is for the unconstrained system. We supply `r1` as
an ArrayPartition of the component functions (in the same order). The arguments of
each component function of `r1` should only take in and return their own parts
of this state.

In the future we can expand this to take in multiple (coupled) systems of this form,
as in FSI problems. In these cases, the component functions should be supplied as
ArrayPartition's of the functions for each system. Each of the component functions of each system
should be able to take in the *full* state and/or constraint and return the
*full* state/constraint, as appropriate, in order to enable coupling of the
systems.
=#


"""
        ConstrainedODEFunction(r1,r2,B1,B2[,L])

This specifies the functions and operators that comprise an ODE problem with the form

``
M\\dfrac{dy}{dt} = Ly - B_1 z + r_1(y,t)
``

``
B_2 y = r_2(t)
``

The optional linear operator `L` defaults to zeros. The `B1` and `B2` functions must be of the respective
in-place forms `B1(dy,z,p)` (to compute the action of `B1` on `z`) and `B2(dz,y,p)` (to compute that action
of `B2` on `y`). The function `r1` must of the in-place form `r1(dy,y,p,t)`, and `r2` must be in the in-place form
`r2(dz,p,t)`.

The optional keyword argument `mass_matrix` can be used to set the mass matrix `M`.

An optional keyword argument `param_update_func` can be used to set a function that updates problem parameters with
the current solution. This function must take the in-place form `f(q,u,p,t)` to update some `q`. (`q` might
simply be `p`.) This function can be used to update `B1` and `B2` with state information, for example.
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

function ConstrainedODEFunction(r1,r2,B1,B2,L=DiffEqLinearOperator(0*I);
                                                 param_update_func = DEFAULT_UPDATE_FUNC,
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



    allinplace(r1,r2,B1,B2) || alloutofplace(r1,r2,B1,B2) || error("Inconsistent function signatures")

    L_local = (L isa DiffEqLinearOperator) ? L : DiffEqLinearOperator(L)

    fill!(_func_cache,0.0)
    odef_nl = SplitFunction(_complete_r1(r1,_func_cache=_func_cache),
                            _complete_B1(B1);_func_cache=deepcopy(_func_cache))
    odef = SplitFunction(L_local, odef_nl ;_func_cache=deepcopy(_func_cache))
    conf = SplitFunction(_complete_r2(r2),_complete_B2(B2);_func_cache=deepcopy(_func_cache))

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

@inline _ode_L!(du,f::ConstrainedODEFunction,u,p,t) = f.odef.f1(du,u,p,t) # L
@inline _ode_r1!(du,f::ConstrainedODEFunction,u,p,t) = f.odef.f2.f.f1(du,u,p,t) # r_1

@inline _ode_neg_B1!(du,f::ConstrainedODEFunction,u,p,t) = f.odef.f2.f.f2(du,u,p,t) # -B_1^T
@inline _constraint_neg_B2!(du,f::ConstrainedODEFunction,u,p,t) = f.conf.f2(du,u,p,t) # -B_2
@inline _constraint_r2!(du,f::ConstrainedODEFunction,u,p,t) = f.conf.f1(du,u,p,t) # r_2

allinplace(r1,r2,B1,B2) = _isinplace_r1(r1) && _isinplace_r2(r2) && _isinplace_B1(B1) && _isinplace_B2(B2)
alloutofplace(r1,r2,B1,B2) = _isoop_r1(r1) && _isoop_r2(r2) && _isoop_B1(B1) && _isoop_B2(B2)

hasaux(r1) = false
hasaux(r1::ArrayPartition) = true

_state_r1(r1) = r1
_state_r1(r1::ArrayPartition) = r1.x[1]
_aux_r1(r1::ArrayPartition) = r1.x[2]


for (f,nv) in ((:r1,4),(:r2,3),(:B1,3),(:B2,3))
  iipfcn = Symbol("_isinplace_",string(f))
  oopfcn = Symbol("_isoop_",string(f))
  completefcn = Symbol("_complete_",string(f))
  @eval $iipfcn($f) = isinplace($f,$nv)
  @eval $iipfcn($f::ArrayPartition) = all(($iipfcn(x) for x in $f.x))
  @eval $oopfcn($f) = numargs($f) == $(nv-1)
  @eval $oopfcn($f::ArrayPartition) = all(($oopfcn(x) for x in $f.x))
  @eval $completefcn($f;_func_cache=nothing) = $completefcn($f,Val($iipfcn($f)),_func_cache)
end

_complete_r1(r1,::Val{true},_func_cache) = (du,u,p,t) -> (dy = state(du); r1(dy,state(u),p,t))
_complete_r1(r1,::Val{false},_func_cache) = (u,p,t) -> r1(state(u),p,t)
_complete_r1(r1::ArrayPartition,::Val{true},_func_cache) =
            SplitFunction((du,u,p,t) ->(dy = state(du); _state_r1(r1)(dy,state(u),p,t)),
                          (du,u,p,t) ->(dy = aux_state(du); _aux_r1(r1)(dy,aux_state(u),p,t));
                          _func_cache=deepcopy(_func_cache))
_complete_r1(r1::ArrayPartition,::Val{false},_func_cache) =
            SplitFunction((du,u,p) ->r1.x[1](state(u),p,t),
                          (du,u,p) ->r1.x[2](aux_state(u),p,t))


_complete_r2(r2,::Val{true},_func_cache) = (du,u,p,t) -> r2(constraint(du),p,t)
_complete_r2(r2,::Val{false},_func_cache) = (u,p,t) -> r2(p,t)


_complete_B1(B1,::Val{true},_func_cache) = (du,u,p,t) -> (dy = state(du); B1(dy,constraint(u),p); dy .*= -1.0)
_complete_B1(B1,::Val{false},_func_cache) = (u,p,t) -> -B1(constraint(u),p)

_complete_B2(B2,::Val{true},_func_cache) = (du,u,p,t) -> (dz = constraint(du); B2(dz,state(u),p); dz .*= -1.0)
_complete_B2(B2,::Val{false},_func_cache) = (u,p,t) -> -B2(state(u),p)

function (f::ConstrainedODEFunction)(du,u,p,t)
    fill!(f.cache,0.0)
    f.odef(f.cache,u,p,t)
    f.conf(du,u,p,t)
    du .+= f.cache
end

@inline isstatic(f::ConstrainedODEFunction) = f.param_update_func == DEFAULT_UPDATE_FUNC ? true : false
