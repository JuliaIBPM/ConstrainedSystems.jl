## Definitions of algorithms ##

# The sc parameter specifies whether it contains static constraint operators or not
# If false, then it expects that the state vector contains a component for updating the opertors

# WrayHERK is scheme C in Liska and Colonius (JCP 2016)
# BraseyHairerHERK is scheme B in Liska and Colonius (JCP 2016)
# LiskaIFHERK is scheme A in Liska and Colonius (JCP 2016)

abstract type ConstrainedOrdinaryDiffEqAlgorithm{sc} <: OrdinaryDiffEq.OrdinaryDiffEqAlgorithm end

for (Alg,Order) in [(:WrayHERK,3),(:BraseyHairerHERK,3),(:LiskaIFHERK,2),(:IFHEEuler,1)]
    @eval struct $Alg{sc,solverType} <: ConstrainedOrdinaryDiffEqAlgorithm{sc}
    end

    @eval $Alg(;static_constraints=true,saddlesolver=Direct) = $Alg{static_constraints,saddlesolver}()

    @eval export $Alg

    @eval alg_order(alg::$Alg) = $Order
end

### Caches ###

# LiskaIFHERK

@cache struct LiskaIFHERKCache{sc,solverType,uType,rateType,expType1,expType2,saddleType,pType,TabType} <: OrdinaryDiffEqMutableCache
  u::uType
  uprev::uType # qi
  k1::rateType # w1
  k2::rateType # w2
  k3::rateType # w3
  utmp::uType  # cache
  dutmp::rateType # cache for rates
  fsalfirst::rateType
  Hhalfdt::expType1
  Hzero::expType2
  S::saddleType
  pnew::pType
  pold::pType
  k::rateType
  tab::TabType
end

struct LiskaIFHERKConstantCache{T,T2} <: OrdinaryDiffEqConstantCache
  ã11::T
  ã21::T
  ã22::T
  ã31::T
  ã32::T
  ã33::T
  c̃1::T2
  c̃2::T2
  c̃3::T2

  function LiskaIFHERKConstantCache(T, T2)
    ã11 = T(1//2)
    ã21 = T(√3/3)
    ã22 = T((3-√3)/3)
    ã31 = T((3+√3)/6)
    ã32 = T(-√3/3)
    ã33 = T((3+√3)/6)
    c̃1 = T2(1//2)
    c̃2 = T2(1.0)
    c̃3 = T2(1.0)
    new{T,T2}(ã11,ã21,ã22,ã31,ã32,ã33,c̃1,c̃2,c̃3)
  end
end

LiskaIFHERKCache{sc,solverType}(u,uprev,k1,k2,k3,utmp,dutmp,fsalfirst,
                                Hhalfdt,Hzero,S,pnew,pold,k,tab) where {sc,solverType} =
        LiskaIFHERKCache{sc,solverType,typeof(u),typeof(k1),typeof(Hhalfdt),typeof(Hzero),
                        typeof(S),typeof(pnew),typeof(tab)}(u,uprev,k1,k2,k3,utmp,dutmp,fsalfirst,
                                                          Hhalfdt,Hzero,S,pnew,pold,k,tab)


function alg_cache(alg::LiskaIFHERK{sc,solverType},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {sc,solverType}

  typeof(u) <: ArrayPartition || error("u must be of type ArrayPartition")

  y, z = u.x[1], u.x[2]

  utmp = zero(u)
  k1, k2, k3, dutmp, fsalfirst, k = (zero(rate_prototype) for i in 1:6)

  tab = LiskaIFHERKConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))

  A = f.odef.f1.f
  Hhalfdt = exp(A,-dt/2,y)
  Hzero = exp(A,zero(dt),y)

  S1 = saddle_system(Hhalfdt,f,p,p,dutmp,solverType)
  S = [S1]
  push!(S,saddle_system(Hzero,f,p,p,dutmp,solverType))

  LiskaIFHERKCache{sc,solverType}(u,uprev,k1,k2,k3,utmp,dutmp,fsalfirst,Hhalfdt,Hzero,S,deepcopy(p),deepcopy(p),k,tab)
end

function alg_cache(alg::LiskaIFHERK,u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{false})
  LiskaIFHERKConstantCache(constvalue(uBottomEltypeNoUnits), constvalue(tTypeNoUnits))
end

function saddle_system(A,f,p,pold,ducache,solver)
    nully, nullz = state(ducache), constraint(ducache)
    B₁ᵀ(z) = (fill!(ducache,0.0); fill!(nully,0.0); -_ode_neg_B1!(ducache,f,ArrayPartition(nully,z),pold,0.0))
    B₂(y) = (fill!(ducache,0.0); fill!(nullz,0.0); -_constraint_neg_B2!(ducache,f,ArrayPartition(y,nullz),p,0.0))
    SaddleSystem(A,B₂,B₁ᵀ,ducache,solver=solver)
end

@inline saddle_system(Sold::SaddleSystem,A,f,p,pold,ducache,solver,::Val{false}) = saddle_system(A,f,p,pold,ducache,solver)
@inline saddle_system(Sold::SaddleSystem,A,f,p,pold,ducache,solver,::Val{true}) = Sold

@inline saddle_system!(S::SaddleSystem,A,f,p,pold,ducache,solver,::Val{false}) = (S = saddle_system(A,f,p,pold,ducache,solver); nothing)
@inline saddle_system!(S::SaddleSystem,A,f,p,pold,ducache,solver,::Val{true}) = nothing

function initialize!(integrator,cache::LiskaIFHERKCache)
    @unpack k,fsalfirst = cache

    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.odef(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    integrator.f.param_update_func(integrator.p,integrator.uprev,integrator.p,integrator.t)
    integrator.destats.nf += 1

end

@muladd function perform_step!(integrator,cache::LiskaIFHERKCache{sc,solverType},repeat_step=false) where {sc,solverType}
    @unpack t,dt,uprev,u,f,p = integrator
    @unpack k1,k2,k3,utmp,dutmp,fsalfirst,Hhalfdt,Hzero,S,pnew,pold,k = cache
    @unpack ã11,ã21,ã22,ã31,ã32,ã33,c̃1,c̃2,c̃3 = cache.tab
    @unpack param_update_func = f

    recursivecopy!(pnew,p)

    # aliases to the state and constraint parts
    ytmp, ztmp = state(utmp), constraint(utmp)
    yprev = state(uprev)
    z = constraint(u)

    ttmp = t
    u .= uprev

    _ode_r1!(k1,f,u,pnew,ttmp)  # should be able to replace this with fsalfirst
    integrator.destats.nf += 1
    @.. k1 *= dt*ã11
    @.. utmp = uprev + k1

    # if applicable, update p, construct new saddle system here, using Hhalfdt
    recursivecopy!(pold,pnew)
    param_update_func(pnew,u,pold,ttmp)
    S[1] = saddle_system(S[1],Hhalfdt,f,pnew,pnew,dutmp,solverType,Val(sc))

    ttmp = t + dt*c̃1
    _constraint_r2!(utmp,f,u,pnew,ttmp) # this should only update the z part
    u .= S[1]\utmp

    ytmp .= typeof(ytmp)(S[1].A⁻¹B₁ᵀf)

    ldiv!(yprev,Hhalfdt,yprev)
    ldiv!(k1.x[1],Hhalfdt,k1.x[1])


    @.. k1 = (k1-utmp)/(dt*ã11)

    _ode_r1!(k2,f,u,pnew,ttmp)
    integrator.destats.nf += 1
    @.. k2 *= dt*ã22
    @.. utmp = uprev + k2 + dt*ã21*k1

    # if applicable, update p, construct new saddle system here, using Hhalfdt
    recursivecopy!(pold,pnew)
    param_update_func(pnew,u,pold,ttmp)
    S[1] = saddle_system(S[1],Hhalfdt,f,pnew,pnew,dutmp,solverType,Val(sc))

    ttmp = t + dt*c̃2
    _constraint_r2!(utmp,f,u,pnew,ttmp)
    u .= S[1]\utmp
    ytmp .= typeof(ytmp)(S[1].A⁻¹B₁ᵀf)

    ldiv!(yprev,Hhalfdt,yprev)
    ldiv!(k1.x[1],Hhalfdt,k1.x[1])
    ldiv!(k2.x[1],Hhalfdt,k2.x[1])

    @.. k2 = (k2-utmp)/(dt*ã22)
    _ode_r1!(k3,f,u,pnew,ttmp)
    integrator.destats.nf += 1
    @.. k3 *= dt*ã33
    @.. utmp = uprev + k3 + dt*ã32*k2 + dt*ã31*k1

    # if applicable, update p, construct new saddle system here, using Hzero (identity)
    recursivecopy!(pold,pnew)
    param_update_func(pnew,u,pold,ttmp)
    S[2] = saddle_system(S[2],Hzero,f,pnew,pnew,dutmp,solverType,Val(sc))

    _constraint_r2!(utmp,f,u,pnew,t+dt)
    u .= S[2]\utmp

    @.. z /= dt*ã33

    param_update_func(pnew,u,p,t)
    f.odef(integrator.fsallast, u, pnew, t+dt)

    recursivecopy!(p,pnew)

    integrator.destats.nf += 1
    return nothing
end
