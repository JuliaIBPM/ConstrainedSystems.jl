## Definitions of algorithms ##


# WrayHERK is scheme C in Liska and Colonius (JCP 2016)
# BraseyHairerHERK is scheme B in Liska and Colonius (JCP 2016)
# LiskaIFHERK is scheme A in Liska and Colonius (JCP 2016)


for (Alg,Order) in [(:WrayHERK,3),(:BraseyHairerHERK,3),(:LiskaIFHERK,2),(:IFHEEuler,1)]
    @eval struct $Alg{solverType} <: ConstrainedOrdinaryDiffEqAlgorithm
      maxiter :: Int
      tol :: Float64
    end

    @eval $Alg(;saddlesolver=Direct,maxiter=4,tol=eps(Float64)) = $Alg{saddlesolver}(maxiter,tol)

    @eval export $Alg

    @eval alg_order(alg::$Alg) = $Order
end



# LiskaIFHERK

@cache struct LiskaIFHERKCache{sc,ni,solverType,uType,rateType,expType1,expType2,saddleType,pType,TabType} <: ConstrainedODEMutableCache{sc,solverType}
  u::uType
  uprev::uType # qi
  k1::rateType # w1
  k2::rateType # w2
  k3::rateType # w3
  utmp::uType  # cache
  udiff::uType
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

struct LiskaIFHERKConstantCache{sc,ni,solverType,T,T2} <: ConstrainedODEConstantCache{sc,solverType}
  ã11::T
  ã21::T
  ã22::T
  ã31::T
  ã32::T
  ã33::T
  c̃1::T2
  c̃2::T2
  c̃3::T2

  function LiskaIFHERKConstantCache{sc,ni,solverType}(T, T2) where {sc,ni,solverType}
    ã11 = T(1//2)
    ã21 = T(√3/3)
    ã22 = T((3-√3)/3)
    ã31 = T((3+√3)/6)
    ã32 = T(-√3/3)
    ã33 = T((3+√3)/6)
    c̃1 = T2(1//2)
    c̃2 = T2(1.0)
    c̃3 = T2(1.0)
    new{sc,ni,solverType,T,T2}(ã11,ã21,ã22,ã31,ã32,ã33,c̃1,c̃2,c̃3)
  end
end

LiskaIFHERKCache{sc,ni,solverType}(u,uprev,k1,k2,k3,utmp,udiff,dutmp,fsalfirst,
                                Hhalfdt,Hzero,S,pnew,pold,k,tab) where {sc,ni,solverType} =
        LiskaIFHERKCache{sc,ni,solverType,typeof(u),typeof(k1),typeof(Hhalfdt),typeof(Hzero),
                        typeof(S),typeof(pnew),typeof(tab)}(u,uprev,k1,k2,k3,utmp,udiff,dutmp,fsalfirst,
                                                          Hhalfdt,Hzero,S,pnew,pold,k,tab)

function alg_cache(alg::LiskaIFHERK{solverType},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {solverType}

  typeof(u) <: ArrayPartition || error("u must be of type ArrayPartition")

  y, z = state(u), constraint(u)

  utmp, udiff = (zero(u) for i in 1:2)
  k1, k2, k3, dutmp, fsalfirst, k = (zero(rate_prototype) for i in 1:6)

  sc = isstatic(f)
  ni = needs_iteration(f,u,p,rate_prototype)

  tab = LiskaIFHERKConstantCache{sc,ni,solverType}(constvalue(uBottomEltypeNoUnits),
                                                constvalue(tTypeNoUnits))

  L = _fetch_ode_L(f)
  Hhalfdt = exp(L,-dt/2,y)
  Hzero = exp(L,zero(dt),y)

  S = []
  push!(S,SaddleSystem(Hhalfdt,f,p,p,dutmp,solverType))
  push!(S,SaddleSystem(Hzero,f,p,p,dutmp,solverType))


  LiskaIFHERKCache{sc,ni,solverType}(u,uprev,k1,k2,k3,utmp,udiff,dutmp,fsalfirst,
                                  Hhalfdt,Hzero,S,deepcopy(p),deepcopy(p),k,tab)
end

function alg_cache(alg::LiskaIFHERK{solverType},u,rate_prototype,
                                  uEltypeNoUnits,uBottomEltypeNoUnits,
                                  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,
                                  p,calck,::Val{false}) where {solverType}
  LiskaIFHERKConstantCache{isstatic(f),
                           needs_iteration(f,u,p,rate_prototype),
                           solverType}(constvalue(uBottomEltypeNoUnits),
                                          constvalue(tTypeNoUnits))
end


# IFHEEuler

@cache struct IFHEEulerCache{sc,ni,solverType,uType,rateType,expType,saddleType,pType} <: ConstrainedODEMutableCache{sc,solverType}
  u::uType
  uprev::uType # qi
  k1::rateType # w1
  utmp::uType  # cache
  udiff::uType
  dutmp::rateType # cache for rates
  fsalfirst::rateType
  Hdt::expType
  S::saddleType
  pnew::pType
  pold::pType
  k::rateType
end

struct IFHEEulerConstantCache{sc,ni,solverType} <: ConstrainedODEConstantCache{sc,solverType}

end

IFHEEulerCache{sc,ni,solverType}(u,uprev,k1,utmp,udiff,dutmp,fsalfirst,
                                Hdt,S,pnew,pold,k) where {sc,ni,solverType} =
        IFHEEulerCache{sc,ni,solverType,typeof(u),typeof(k1),typeof(Hdt),
                        typeof(S),typeof(pnew)}(u,uprev,k1,utmp,udiff,dutmp,fsalfirst,
                                                              Hdt,S,pnew,pold,k)

function alg_cache(alg::IFHEEuler{solverType},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {solverType}

  typeof(u) <: ArrayPartition || error("u must be of type ArrayPartition")

  y, z = state(u), constraint(u)

  utmp, udiff = (zero(u) for i in 1:2)
  k1, dutmp, fsalfirst, k = (zero(rate_prototype) for i in 1:4)

  sc = isstatic(f)
  ni = needs_iteration(f,u,p,rate_prototype)

  Hdt = exp(_fetch_ode_L(f),-dt,y)

  S = []
  push!(S,SaddleSystem(Hdt,f,p,p,dutmp,solverType))

  IFHEEulerCache{sc,ni,solverType}(u,uprev,k1,utmp,udiff,dutmp,fsalfirst,
                                  Hdt,S,deepcopy(p),deepcopy(p),k)
end

function alg_cache(alg::IFHEEuler{solverType},u,rate_prototype,
                                  uEltypeNoUnits,uBottomEltypeNoUnits,
                                  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,
                                  p,calck,::Val{false}) where {solverType}
  IFHEEulerConstantCache{isstatic(f),needs_iteration(f,u,p,rate_prototype),solverType}()
end


#######

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


@muladd function perform_step!(integrator,cache::LiskaIFHERKCache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}
    @unpack t,dt,uprev,u,f,p,alg = integrator
    @unpack maxiter, tol = alg
    @unpack k1,k2,k3,utmp,udiff,dutmp,fsalfirst,Hhalfdt,Hzero,S,pnew,pold,k = cache
    @unpack ã11,ã21,ã22,ã31,ã32,ã33,c̃1,c̃2,c̃3 = cache.tab
    @unpack param_update_func = f

    init_err = float(1)
    init_iter = ni ? 1 : maxiter

    # aliases to the state and constraint parts
    ytmp, ztmp, xtmp = state(utmp), constraint(utmp), aux_state(utmp)
    yprev = state(uprev)
    y, z, x = state(u), constraint(u), aux_state(u)

    ttmp = t
    u .= uprev
    recursivecopy!(pnew,p)

    ## Stage 1
    _ode_r1!(k1,f,u,pnew,ttmp)
    integrator.destats.nf += 1
    @.. k1 *= dt*ã11
    @.. utmp = uprev + k1
    ttmp = t + dt*c̃1

    # if applicable, update p, construct new saddle system here, using Hhalfdt
    # and solve system. Solve iteratively if saddle operators depend on
    # constrained part of the state.
    recursivecopy!(pold,pnew)
    err, numiter = init_err, init_iter
    u .= utmp # initial guess for iterations
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew,u,pold,ttmp)
      S[1] = SaddleSystem(S[1],Hhalfdt,f,pnew,pold,cache)
      _constraint_r2!(utmp,f,u,pnew,ttmp) # only updates the z part
      mainvector(u) .= S[1]\mainvector(utmp)
      @.. udiff -= u
      numiter += 1
      err = _l2norm(udiff)
    end
    zero_vec!(xtmp)
    B1_times_z!(utmp,S[1])

    ldiv!(yprev,Hhalfdt,yprev)
    ldiv!(state(k1),Hhalfdt,state(k1))

    @.. k1 = (k1-utmp)/(dt*ã11)  # r1(y,t) - B1T*z

    ## Stage 2
    _ode_r1!(k2,f,u,pnew,ttmp)
    integrator.destats.nf += 1
    @.. k2 *= dt*ã22
    @.. utmp = uprev + k2 + dt*ã21*k1
    ttmp = t + dt*c̃2

    # if applicable, update p, construct new saddle system here, using Hhalfdt
    recursivecopy!(pold,pnew)
    err, numiter = init_err, init_iter
    u .= utmp
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew,u,pold,ttmp)
      S[1] = SaddleSystem(S[1],Hhalfdt,f,pnew,pold,cache)
      _constraint_r2!(utmp,f,u,pnew,ttmp)
      mainvector(u) .= S[1]\mainvector(utmp)
      @.. udiff -= u
      numiter += 1
      err = _l2norm(udiff)
    end
    zero_vec!(xtmp)
    B1_times_z!(utmp,S[1])

    ldiv!(yprev,Hhalfdt,yprev)
    ldiv!(state(k1),Hhalfdt,state(k1))
    ldiv!(state(k2),Hhalfdt,state(k2))

    @.. k2 = (k2-utmp)/(dt*ã22)

    ## Stage 3
    _ode_r1!(k3,f,u,pnew,ttmp)
    integrator.destats.nf += 1
    @.. k3 *= dt*ã33
    @.. utmp = uprev + k3 + dt*ã32*k2 + dt*ã31*k1
    ttmp = t + dt

    # if applicable, update p, construct new saddle system here, using Hzero (identity)
    recursivecopy!(pold,pnew)
    err, numiter = init_err, init_iter
    u .= utmp
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew,u,pold,ttmp)
      S[2] = SaddleSystem(S[2],Hzero,f,pnew,pold,cache)
      _constraint_r2!(utmp,f,u,pnew,t+dt)
      mainvector(u) .= S[2]\mainvector(utmp)
      @.. udiff -= u
      numiter += 1
      err = _l2norm(udiff)
      #println("error = ",err)
    end

    @.. z /= (dt*ã33)

    param_update_func(pnew,u,pold,t+dt)

    f.odef(integrator.fsallast, u, pnew, t+dt)

    recursivecopy!(p,pnew)

    integrator.destats.nf += 1
    return nothing
end

####

function initialize!(integrator,cache::IFHEEulerCache)
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

@muladd function perform_step!(integrator,cache::IFHEEulerCache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}
    @unpack t,dt,uprev,u,f,p,alg = integrator
    @unpack k1,utmp,udiff,dutmp,fsalfirst,Hdt,S,pnew,pold,k = cache
    @unpack maxiter, tol = alg
    @unpack param_update_func = f

    init_err = float(1)
    init_iter = ni ? 1 : maxiter

    # aliases to the state and constraint parts
    ytmp, ztmp = state(utmp), constraint(utmp)
    z = constraint(u)

    ttmp = t
    u .= uprev
    recursivecopy!(pnew,p)

    _ode_r1!(k1,f,u,pnew,ttmp)
    integrator.destats.nf += 1
    @.. k1 *= dt
    @.. utmp = uprev + k1
    ttmp = t + dt

    # if applicable, update p, construct new saddle system here, using Hdt
    recursivecopy!(pold,pnew)
    err, numiter = init_err, init_iter
    u .= utmp
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew,u,pold,ttmp)

      S[1] = SaddleSystem(S[1],Hdt,f,pnew,pold,cache)

      _constraint_r2!(utmp,f,u,pnew,t+dt)

      mainvector(u) .= S[1]\mainvector(utmp)
      @.. udiff -= u
      numiter += 1
      err = _l2norm(udiff)
      #println("error = ",err)
    end

    @.. z /= dt

    param_update_func(pnew,u,p,t)
    f.odef(integrator.fsallast, u, pnew, t+dt)

    recursivecopy!(p,pnew)

    integrator.destats.nf += 1
    return nothing
end


function initialize!(integrator,cache::IFHEEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f.odef(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.p = integrator.f.param_update_func(integrator.uprev,integrator.p,integrator.t)
  integrator.destats.nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::IFHEEulerConstantCache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}
  @unpack t,dt,uprev,f,p, alg = integrator
  @unpack maxiter, tol = alg
  @unpack param_update_func = f

  udiff = deepcopy(uprev)
  ducache = deepcopy(uprev)

  init_err = float(1)
  init_iter = ni ? 1 : maxiter

  L = _fetch_ode_L(f)
  Hdt = exp(L,-dt,state(uprev))

  pnew = deepcopy(p)

  k1 = _ode_r1(f,uprev,pnew,t)
  integrator.destats.nf += 1
  @.. k1 *= dt
  utmp = @.. uprev + k1

  # if applicable, update p, construct new saddle system here, using Hdt
  pold = deepcopy(pnew)
  err, numiter = init_err, init_iter
  u = deepcopy(utmp)
  while err > tol && numiter <= maxiter
    udiff .= u
    pnew = param_update_func(u,pold,t+dt)
    S = SaddleSystem(Hdt,f,pnew,pold,ducache,solverType)
    constraint(utmp) .= constraint(_constraint_r2(f,u,pnew,t+dt))
    mainvector(u) .= S\mainvector(utmp)
    @.. udiff -= u
    numiter += 1
    err = _l2norm(udiff)
    #println("error = ",err)
  end

  @.. constraint(u) /= dt

  pnew = param_update_func(u,p,t)
  k = f.odef(u, pnew, t+dt)
  recursivecopy!(p,pnew)

  integrator.destats.nf += 1
  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
