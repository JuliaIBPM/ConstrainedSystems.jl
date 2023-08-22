export DiffEqLinearOperator, ConstrainedODEFunction


DEFAULT_PARAM_UPDATE_FUNC(q,u,p,t) = q
DEFAULT_PARAM_UPDATE_FUNC(u,p,t) = p


### Abstract algorithms ###
abstract type ConstrainedOrdinaryDiffEqAlgorithm <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

### Abstract caches ###
# The sc parameter specifies whether it contains static constraint operators or not
# If false, then it expects that the state vector contains a component for updating the opertors
abstract type ConstrainedODEMutableCache{sc,solverType} <: OrdinaryDiffEqMutableCache end
abstract type ConstrainedODEConstantCache{sc,solverType} <: OrdinaryDiffEqConstantCache end


#### Operator and function types ####

mutable struct DiffEqLinearOperator{T,aType} <: AbstractDiffEqLinearOperator{T}
    L :: aType
    DiffEqLinearOperator(L::aType; update_func=DEFAULT_PARAM_UPDATE_FUNC,
                          dtype=Float64) where {aType} = new{dtype,aType}(L)
end

(f::DiffEqLinearOperator)(du,u,p,t) = (dy = state(du); mul!(dy,f.L,state(u)))
(f::DiffEqLinearOperator)(u,p,t) = (du = deepcopy(u); zero_vec!(du); dy = state(du); mul!(dy,f.L,state(u)); return du)


import Base: exp
exp(f::DiffEqLinearOperator,args...) = exp(f.L,args...)
exp(f::ArrayPartition,t,x::ArrayPartition) =
                ArrayPartition((exp(Li,t,xi) for (Li,xi) in zip(f.L.x,x.x))...)

exp(f::ODEFunction,args...) = exp(f.f,args...)


has_exp(::DiffEqLinearOperator) = true

exp(L::AbstractMatrix,t,x) = exp(factorize(L)*t)
exp(L::UniformScaling,t,x) = exp(Diagonal(L,length(x))*t)

#=
ConstrainedODEFunction

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
        ConstrainedODEFunction(r1,r2,B1,B2[,L][,C])

This specifies the functions and operators that comprise an ODE problem with the form

``
\\dfrac{dy}{dt} = Ly - B_1 z + r_1(y,t)
``

``
B_2 y + C z = r_2(x,t)
``

where ``y`` is the state, ``z`` is a constraint force, and ``x`` is an auxiliary
state describing the constraints.

The optional linear operator `L` defaults to zeros. The `B1` and `B2` functions must be of the respective
in-place forms `B1(dy,z,x,p)` (to compute the action of `B1` on `z`) and `B2(dz,y,x,p)` (to compute the action
of `B2` on `y`). The function `r1` must of the in-place form `r1(dy,y,x,p,t)`, and `r2` must be in the in-place form
`r2(dz,x,p,t)`. The `C` function can be omitted, but if it is included, then it must be of the form
`C(dz,z,x,p)` (to compute the action of `C` on `z`). Alternatively, one can supply out-of-place forms, respectively, as `B1(z,x,p)`, `B2(y,x,p)`,
`C(z,x,p)`, `r1(y,x,p,t)` and `r2(x,p,t)`.

An optional keyword argument `param_update_func` can be used to set a function that updates problem parameters with
the current solution. This function must take the in-place form `f(q,u,p,t)` or out of place form
`f(u,p,t)` to create some `q` based on `u`, where `y = state(u)`, `z = constraint(u)` and `x = aux_state(u)`.
(Note that `q` might enter the function simply as `p`, to be mutated.) This function can
be used to update `B1`, `B2`, and `C`, for example.

We can also include another (unconstrained) set of equations to the set above
in order to update `x`:

``
\\dfrac{dx}{dt} = r_{1x}(u,p,t)
``

In this case, the right-hand side has access to the entire `u` vector. We would pass
the pair of `r1` functions as an `ArrayPartition`.
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

function ConstrainedODEFunction(r1,r2,B1,B2,L=DiffEqLinearOperator(0*I),C=nothing;
                                                 param_update_func = DEFAULT_PARAM_UPDATE_FUNC,
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



    allempty(r2,B1,B2) || noneempty(r2,B1,B2) || error("Inconsistent null operators")
    unconstrained = allempty(r2,B1,B2)

    allinplace(r1,r2,B1,B2) || alloutofplace(r1,r2,B1,B2) || error("Inconsistent function signatures")
    iip = allinplace(r1,r2,B1,B2)

    static = param_update_func == DEFAULT_PARAM_UPDATE_FUNC ? true : false


    L_local = (L isa DiffEqLinearOperator) ? L : DiffEqLinearOperator(L)

    local_cache = deepcopy(_func_cache)
    zero_vec!(local_cache)

    odef_nl = SplitFunction(_complete_r1(r1,Val(iip),_func_cache=local_cache),
                            _complete_B1(B1,Val(iip));_func_cache=deepcopy(local_cache))
    odef = SplitFunction(L_local, odef_nl ;_func_cache=deepcopy(local_cache))
    conf_lhs = SplitFunction(_complete_B2(B2,Val(iip)),_complete_C(C,Val(iip));_func_cache=deepcopy(local_cache))
    conf = SplitFunction(_complete_r2(r2,Val(iip)),conf_lhs;_func_cache=deepcopy(local_cache))

    ConstrainedODEFunction{iip,static,typeof(odef),
                           typeof(conf),typeof(mass_matrix),typeof(local_cache),
                           typeof(analytic),typeof(tgrad),typeof(jac),typeof(jvp),typeof(vjp),
                           typeof(jac_prototype),typeof(sparsity),
                           typeof(Wfact),typeof(Wfact_t),typeof(paramjac),typeof(syms),
                           typeof(colorvec),typeof(param_update_func)}(odef,conf,mass_matrix,local_cache,
                           analytic,tgrad,jac,jvp,vjp,jac_prototype,
                           sparsity,Wfact,Wfact_t,paramjac,syms,colorvec,param_update_func)
end

# For unconstrained systems
ConstrainedODEFunction(r1,L=DiffEqLinearOperator(0*I);kwargs...) =
      ConstrainedODEFunction(r1,nothing,nothing,nothing,L,nothing;kwargs...)


function Base.show(io::IO, m::MIME"text/plain",f::ConstrainedODEFunction{iip,static}) where {iip,static}
    iips = iip ? "in-place" : "out-of-place"
    statics = static ? "static" : "variable"
    println(io,"Constrained ODE function of $iips type and $statics constraints")
end

# Here is where we define the structure of the function
@inline _fetch_ode_L(f::ConstrainedODEFunction) = f.odef.f1
@inline _fetch_ode_r1(f::ConstrainedODEFunction) = f.odef.f2.f.f1
@inline _fetch_ode_neg_B1(f::ConstrainedODEFunction) = f.odef.f2.f.f2
@inline _fetch_constraint_r2(f::ConstrainedODEFunction) = f.conf.f1
@inline _fetch_constraint_neg_B2(f::ConstrainedODEFunction) = f.conf.f2.f.f1
@inline _fetch_constraint_neg_C(f::ConstrainedODEFunction) = f.conf.f2.f.f2


for fcn in (:_ode_L,:_ode_r1,:_ode_neg_B1,:_constraint_neg_B2,:_constraint_neg_C,:_constraint_r2)
  fetchfcn = Symbol("_fetch",string(fcn))
  iipfcn = Symbol(string(fcn),"!")
  @eval $iipfcn(du,f::ConstrainedODEFunction,u,p,t) = $fetchfcn(f)(du,u,p,t)
  @eval $fcn(f::ConstrainedODEFunction,u,p,t) = $fetchfcn(f)(u,p,t)
end

allempty(::Nothing,::Nothing,::Nothing) = true
allempty(r2,B1,B2) = false
noneempty(r2,B1,B2) = !(isnothing(r2) || isnothing(B1) || isnothing(B2))


allinplace(r1,r2,B1,B2) = _isinplace_r1(r1) && _isinplace_r2(r2) && _isinplace_B1(B1) && _isinplace_B2(B2)
alloutofplace(r1,r2,B1,B2) = _isoop_r1(r1) && _isoop_r2(r2) && _isoop_B1(B1) && _isoop_B2(B2)

allinplace(r1,::Nothing,::Nothing,::Nothing) = _isinplace_r1(r1)
alloutofplace(r1,::Nothing,::Nothing,::Nothing) = _isoop_r1(r1)



for (f,nv,nvaux) in ((:r1,5,4),(:r2,4,0),(:B1,4,0),(:B2,4,0),(:C,4,0))
  iipfcn = Symbol("_isinplace_",string(f))
  oopfcn = Symbol("_isoop_",string(f))
  completefcn = Symbol("_complete_",string(f))
  @eval $iipfcn(fcn) = isinplace(fcn,$nv)
  @eval $iipfcn(fcn::ArrayPartition) = $iipfcn(fcn.x[1]) && isinplace(fcn.x[2],$nvaux)
  @eval $oopfcn(fcn) = first(numargs(fcn)) == $(nv-1)
  @eval $oopfcn(fcn::ArrayPartition) = $oopfcn(fcn.x[1]) && first(numargs(fcn.x[2])) == $(nvaux-1)
  @eval $completefcn(fcn,::Val{iip};_func_cache=nothing) where {iip} = $completefcn(fcn,Val(iip),_func_cache)
  @eval $completefcn(::Nothing,::Val{false},_func_cache) = (u,p,t) -> zero(u)
end

_complete_r1(r1,::Val{true},_func_cache) = (du,u,p,t) -> (dy = state(du); y = state(u); x = aux_state(u); r1(dy,y,x,p,t))
_complete_r1(r1,::Val{false},_func_cache) = (u,p,t) -> (du = deepcopy(u); zero_vec!(du); y = state(u); x = aux_state(u); state(du) .= r1(y,x,p,t); return du)
_complete_r1(r1::ArrayPartition,::Val{true},_func_cache) =
            SplitFunction((du,u,p,t) ->(dy = state(du); dx = aux_state(du); zero_vec!(dx); y = state(u); x = aux_state(u); state_r1(r1)(dy,y,x,p,t)),
                          (du,u,p,t) ->(dy = state(du); dx = aux_state(du); zero_vec!(dy); aux_r1(r1)(dx,u,p,t));
                          _func_cache=deepcopy(_func_cache))
_complete_r1(r1::ArrayPartition,::Val{false},_func_cache) =
            SplitFunction((u,p,t) -> (du = deepcopy(u); zero_vec!(du); y = state(u); x = aux_state(u); state(du) .= state_r1(r1)(y,x,p,t); return du),
                          (u,p,t) -> (du = deepcopy(u); zero_vec!(du); aux_state(du) .= aux_r1(r1)(u,p,t); return du))


_complete_r2(r2,::Val{true},_func_cache) = (du,u,p,t) -> (dz = constraint(du); x = aux_state(u); r2(dz,x,p,t))
_complete_r2(::Nothing,::Val{true},_func_cache) = (du,u,p,t) -> zero_vec!(constraint(du))
_complete_r2(r2,::Val{false},_func_cache) = (u,p,t) -> (du = deepcopy(u); zero_vec!(du); x = aux_state(u); constraint(du) .= r2(x,p,t); return du)


_complete_B1(B1,::Val{true},_func_cache) = (du,u,p,t) -> (dy = state(du); dx = aux_state(du); zero_vec!(dx);
                                              z = constraint(u); x = aux_state(u); B1(dy,z,x,p); dy .*= -1.0)
_complete_B1(::Nothing,::Val{true},_func_cache) = (du,u,p,t) -> (zero_vec!(aux_state(du)); zero_vec!(state(du)))
_complete_B1(B1,::Val{false},_func_cache) = (u,p,t) -> (du = deepcopy(u); zero_vec!(du);
                                              z = constraint(u); x = aux_state(u); state(du) .= -B1(z,x,p); return du)

_complete_B2(B2,::Val{true},_func_cache) = (du,u,p,t) -> (dz = constraint(du); y = state(u); x = aux_state(u);
                                              B2(dz,y,x,p); dz .*= -1.0)
_complete_B2(::Nothing,::Val{true},_func_cache) = (du,u,p,t) -> zero_vec!(constraint(du))
_complete_B2(B2,::Val{false},_func_cache) = (u,p,t) -> (du = deepcopy(u); zero_vec!(du); y = state(u); x = aux_state(u);
                                             constraint(du) .= -B2(y,x,p); return du)


_complete_C(C,::Val{true},_func_cache) = (du,u,p,t) -> (dz = constraint(du); z = constraint(u); x = aux_state(u);
                                              C(dz,z,x,p); dz .*= -1.0)
_complete_C(::Nothing,::Val{true},_func_cache) = (du,u,p,t) -> zero_vec!(constraint(du))
_complete_C(C,::Val{false},_func_cache) = (u,p,t) -> (du = deepcopy(u); zero_vec!(du); z = constraint(u);
                                              x = aux_state(u); constraint(du) .= -C(z,x,p); return du)



function (f::ConstrainedODEFunction)(du,u,p,t)
    zero_vec!(f.cache)
    zero_vec!(du)
    f.odef(f.cache,u,p,t)
    f.conf(du,u,p,t)
    du .+= f.cache
end

(f::ConstrainedODEFunction)(u,p,t) = f.odef(u,p,t) + f.conf(u,p,t)


@inline isstatic(f::ConstrainedODEFunction) = f.param_update_func == DEFAULT_PARAM_UPDATE_FUNC ? true : false
