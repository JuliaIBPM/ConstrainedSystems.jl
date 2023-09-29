## Definitions of algorithms ##

stats_field(integrator) = integrator.stats

# WrayHERK is scheme C in Liska and Colonius (JCP 2016)
# BraseyHairerHERK is scheme B in Liska and Colonius (JCP 2016)
# LiskaIFHERK is scheme A in Liska and Colonius (JCP 2016)
# IFHEEuler is an Euler method with integrating factor
# HETrapezoidalAB2 is a Crank-Nicolson/2nd-order Adams Bashforth method


for (Alg,Order) in [(:WrayHERK,3),(:BraseyHairerHERK,3),(:LiskaIFHERK,2),(:IFHEEuler,1),(:HETrapezoidalAB2,2)]
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
  ptmp::pType
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
                                Hhalfdt,Hzero,S,ptmp,k,tab) where {sc,ni,solverType} =
        LiskaIFHERKCache{sc,ni,solverType,typeof(u),typeof(k1),typeof(Hhalfdt),typeof(Hzero),
                        typeof(S),typeof(ptmp),typeof(tab)}(u,uprev,k1,k2,k3,utmp,udiff,dutmp,fsalfirst,
                                                          Hhalfdt,Hzero,S,ptmp,k,tab)

function alg_cache(alg::LiskaIFHERK{solverType},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {solverType}

  u isa ArrayPartition || error("u must be of type ArrayPartition")

  y, z = state(u), constraint(u)

  utmp, udiff = (zero(u) for i in 1:2)
  k1, k2, k3, dutmp, fsalfirst, k = (zero(rate_prototype) for i in 1:6)

  sc = isstatic(f)
  ni = needs_iteration(f,u,p,rate_prototype)

  tab = LiskaIFHERKConstantCache{sc,ni,solverType}(constvalue(uBottomEltypeNoUnits),
                                                constvalue(tTypeNoUnits))

  @unpack ã11,ã22,ã33 = tab

  L = _fetch_ode_L(f)
  Hhalfdt = exp(L,-dt/2,y)
  Hzero = exp(L,zero(dt),y)

  S = []
  push!(S,SaddleSystem(Hhalfdt,f,p,p,dutmp,solverType;cfact=1.0/(ã11*dt)))
  push!(S,SaddleSystem(Hhalfdt,f,p,p,dutmp,solverType;cfact=1.0/(ã22*dt)))
  push!(S,SaddleSystem(Hzero,f,p,p,dutmp,solverType;cfact=1.0/(ã33*dt)))


  LiskaIFHERKCache{sc,ni,solverType}(u,uprev,k1,k2,k3,utmp,udiff,dutmp,fsalfirst,
                                  Hhalfdt,Hzero,S,deepcopy(p),k,tab)
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
  ptmp::pType
  k::rateType
end

struct IFHEEulerConstantCache{sc,ni,solverType} <: ConstrainedODEConstantCache{sc,solverType}

end

IFHEEulerCache{sc,ni,solverType}(u,uprev,k1,utmp,udiff,dutmp,fsalfirst,
                                Hdt,S,ptmp,k) where {sc,ni,solverType} =
        IFHEEulerCache{sc,ni,solverType,typeof(u),typeof(k1),typeof(Hdt),
                        typeof(S),typeof(ptmp)}(u,uprev,k1,utmp,udiff,dutmp,fsalfirst,
                                                              Hdt,S,ptmp,k)

function alg_cache(alg::IFHEEuler{solverType},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {solverType}

  u isa ArrayPartition || error("u must be of type ArrayPartition")

  y, z = state(u), constraint(u)

  utmp, udiff = (zero(u) for i in 1:2)
  k1, dutmp, fsalfirst, k = (zero(rate_prototype) for i in 1:4)

  sc = isstatic(f)
  ni = needs_iteration(f,u,p,rate_prototype)

  Hdt = exp(_fetch_ode_L(f),-dt,y)

  S = []
  push!(S,SaddleSystem(Hdt,f,p,p,dutmp,solverType;cfact=1.0/dt))

  IFHEEulerCache{sc,ni,solverType}(u,uprev,k1,utmp,udiff,dutmp,fsalfirst,
                                  Hdt,S,deepcopy(p),k)
end

function alg_cache(alg::IFHEEuler{solverType},u,rate_prototype,
                                  uEltypeNoUnits,uBottomEltypeNoUnits,
                                  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,
                                  p,calck,::Val{false}) where {solverType}
  IFHEEulerConstantCache{isstatic(f),needs_iteration(f,u,p,rate_prototype),solverType}()
end

# Half-explicit trapezoidal-Adams/Bashforth 2 (HETrapezoidalAB2)

@cache struct HETrapezoidalAB2Cache{sc,ni,solverType,uType,rateType,implicitType,saddleType,pType,TabType} <: ConstrainedODEMutableCache{sc,solverType}
  u::uType
  uprev::uType
  ki::rateType
  ke::rateType
  utmp::uType  # cache
  udiff::uType
  dutmp::rateType # cache for rates
  fsalfirst::rateType
  A::implicitType
  S::saddleType
  ptmp::pType
  k::rateType
  tab::TabType
end

struct HETrapezoidalAB2ConstantCache{sc,ni,solverType,T} <: ConstrainedODEConstantCache{sc,solverType}
  α̃1::T
  α̃2::T
  β̃1::T
  β̃2::T
  function HETrapezoidalAB2ConstantCache{sc,ni,solverType}(T) where {sc,ni,solverType}
    α̃1 = T(1//2)
    α̃2 = T(1//2)
    β̃1 = T(3//2)
    β̃2 = T(-1//2)
    new{sc,ni,solverType,T}(α̃1,α̃2,β̃1,β̃2)
  end
end

HETrapezoidalAB2Cache{sc,ni,solverType}(u,uprev,ki,ke,utmp,udiff,dutmp,fsalfirst,
                                A,S,ptmp,k,tab) where {sc,ni,solverType} =
        HETrapezoidalAB2Cache{sc,ni,solverType,typeof(u),typeof(ki),typeof(A),
                        typeof(S),typeof(ptmp),typeof(tab)}(u,uprev,ki,ke,utmp,udiff,dutmp,fsalfirst,A,S,ptmp,k,tab)



function alg_cache(alg::HETrapezoidalAB2{solverType},u,rate_prototype,uEltypeNoUnits,uBottomEltypeNoUnits,
                   tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,p,calck,::Val{true}) where {solverType}

  u isa ArrayPartition || error("u must be of type ArrayPartition")

  y, z = state(u), constraint(u)

  utmp, udiff = (zero(u) for i in 1:2)
  ki, ke, dutmp, fsalfirst, k = (zero(rate_prototype) for i in 1:5)

  sc = isstatic(f)
  ni = needs_iteration(f,u,p,rate_prototype)

  tab = HETrapezoidalAB2ConstantCache{sc,ni,solverType}(constvalue(uBottomEltypeNoUnits))

  @unpack α̃1 = tab

  A = implicit_operator(_fetch_ode_L(f),α̃1*dt)

  S = []
  push!(S,SaddleSystem(A,f,p,p,dutmp,solverType;cfact=1.0/(α̃1*dt)))
  push!(S,SaddleSystem(A,f,p,p,dutmp,solverType;cfact=1.0/dt))


  HETrapezoidalAB2Cache{sc,ni,solverType}(u,uprev,ki,ke,utmp,udiff,dutmp,fsalfirst,
                                  A,S,deepcopy(p),k,tab)
end

function alg_cache(alg::HETrapezoidalAB2{solverType},u,rate_prototype,
                                  uEltypeNoUnits,uBottomEltypeNoUnits,
                                  tTypeNoUnits,uprev,uprev2,f,t,dt,reltol,
                                  p,calck,::Val{false}) where {solverType}
  HETrapezoidalAB2ConstantCache{isstatic(f),needs_iteration(f,u,p,rate_prototype),solverType}(constvalue(uBottomEltypeNoUnits))
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
    stats_field(integrator).nf += 1

end


@muladd function perform_step!(integrator,cache::LiskaIFHERKCache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}


    @unpack t,dt,uprev,u,f,p,alg,opts = integrator
    @unpack internalnorm = opts
    @unpack maxiter, tol = alg
    @unpack k1,k2,k3,utmp,udiff,dutmp,fsalfirst,Hhalfdt,Hzero,S,ptmp,k = cache
    @unpack ã11,ã21,ã22,ã31,ã32,ã33,c̃1,c̃2,c̃3 = cache.tab
    @unpack param_update_func = f


    init_err = float(1)
    init_iter = ni ? 1 : maxiter

    # aliases to the state and constraint parts
    ytmp, ztmp, xtmp = state(utmp), constraint(utmp), aux_state(utmp)
    yprev = state(uprev)
    y, z, x = state(u), constraint(u), aux_state(u)
    pold_ptr = p
    pnew_ptr = ptmp

    ttmp = t
    @.. u = uprev


    ## Stage 1
    _ode_full_rhs!(k1,f,u,pold_ptr,ttmp)
    stats_field(integrator).nf += 1
    @.. k1 *= dt*ã11
    @.. utmp = uprev + k1
    ttmp = t + dt*c̃1

    @.. u = utmp

    # if applicable, update p, construct new saddle system here, using Hhalfdt
    # and solve system. Solve iteratively if saddle operators depend on
    # constrained part of the state.
    err, numiter = init_err, init_iter
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew_ptr,u,pold_ptr,ttmp)
      S[1] = SaddleSystem(S[1],Hhalfdt,f,pnew_ptr,pold_ptr,cache)
      _constraint_r2!(utmp,f,u,pnew_ptr,ttmp) # only updates the z part
      ldiv!(mainvector(u),S[1],mainvector(utmp))
      @.. udiff -= u
      numiter += 1
      err = internalnorm(udiff,ttmp)
    end
    zero_vec!(xtmp)
    B1_times_z!(utmp,S[1])
    pold_ptr = ptmp
    pnew_ptr = p

    ldiv!(yprev,Hhalfdt,yprev)
    ldiv!(state(k1),Hhalfdt,state(k1))

    @.. k1 = (k1-utmp)/(dt*ã11)  # r1(y,t) - B1T*z


    ## Stage 2
    _ode_full_rhs!(k2,f,u,pold_ptr,ttmp)
    stats_field(integrator).nf += 1
    @.. k2 *= dt*ã22
    @.. utmp = uprev + k2 + dt*ã21*k1
    ttmp = t + dt*c̃2

    @.. u = utmp

    # if applicable, update p, construct new saddle system here, using Hhalfdt
    err, numiter = init_err, init_iter
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew_ptr,u,pold_ptr,ttmp)
      S[2] = SaddleSystem(S[2],Hhalfdt,f,pnew_ptr,pold_ptr,cache)
      _constraint_r2!(utmp,f,u,pnew_ptr,ttmp)
      ldiv!(mainvector(u),S[2],mainvector(utmp))
      @.. udiff -= u
      numiter += 1
      err = internalnorm(udiff,ttmp)
    end
    zero_vec!(xtmp)
    B1_times_z!(utmp,S[2])
    pold_ptr = p
    pnew_ptr = ptmp

    ldiv!(yprev,Hhalfdt,yprev)
    ldiv!(state(k1),Hhalfdt,state(k1))
    ldiv!(state(k2),Hhalfdt,state(k2))

    @.. k2 = (k2-utmp)/(dt*ã22)


    ## Stage 3
    _ode_full_rhs!(k3,f,u,pold_ptr,ttmp)
    stats_field(integrator).nf += 1
    @.. k3 *= dt*ã33
    @.. utmp = uprev + k3 + dt*ã32*k2 + dt*ã31*k1
    ttmp = t + dt

    @.. u = utmp

    # if applicable, update p, construct new saddle system here, using Hzero (identity)
    err, numiter = init_err, init_iter
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew_ptr,u,pold_ptr,ttmp)
      S[3] = SaddleSystem(S[3],Hzero,f,pnew_ptr,pold_ptr,cache)
      _constraint_r2!(utmp,f,u,pnew_ptr,t+dt)
      ldiv!(mainvector(u),S[3],mainvector(utmp))
      @.. udiff -= u
      numiter += 1
      err = internalnorm(udiff,ttmp)
      #println("error = ",err)
    end
    @.. z /= (dt*ã33)

    # Final steps
    param_update_func(p,u,ptmp,t+dt)
    f.odef(integrator.fsallast, u, p, t+dt)
    stats_field(integrator).nf += 1


    return nothing
end

function initialize!(integrator,cache::LiskaIFHERKConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f.odef(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.p = integrator.f.param_update_func(integrator.uprev,integrator.p,integrator.t)
  stats_field(integrator).nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::LiskaIFHERKConstantCache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}
    @unpack t,dt,uprev,f,p,opts,alg = integrator
    @unpack internalnorm = opts
    @unpack maxiter, tol = alg
    @unpack ã11,ã21,ã22,ã31,ã32,ã33,c̃1,c̃2,c̃3 = cache
    @unpack param_update_func = f

    init_err = float(1)
    init_iter = ni ? 1 : maxiter

    # set up some cache variables
    yprev = state(uprev)
    udiff = deepcopy(uprev)
    ducache = deepcopy(uprev)
    ptmp = deepcopy(p)
    L = _fetch_ode_L(f)
    Hhalfdt = exp(L,-dt/2,state(uprev))
    Hzero = exp(L,zero(dt),state(uprev))
    pold_ptr = p

    ## Stage 1
    ttmp = t
    k1 = _ode_full_rhs(f,uprev,pold_ptr,ttmp)
    stats_field(integrator).nf += 1
    @.. k1 *= dt*ã11
    utmp = @.. uprev + k1
    ttmp = t + dt*c̃1

    # if applicable, update p, construct new saddle system here, using Hhalfdt
    # and solve system. Solve iteratively if saddle operators depend on
    # constrained part of the state.
    err, numiter = init_err, init_iter
    u = deepcopy(utmp)
    while err > tol && numiter <= maxiter
      udiff .= u
      ptmp = param_update_func(u,pold_ptr,ttmp)
      pnew_ptr = ptmp
      S = SaddleSystem(Hhalfdt,f,pnew_ptr,pold_ptr,ducache,solverType;cfact=1.0/(ã11*dt))
      constraint(utmp) .= constraint(_constraint_r2(f,u,pnew_ptr,ttmp))
      ldiv!(mainvector(u),S,mainvector(utmp))
      B1_times_z!(ducache,S)

      @.. udiff -= u
      numiter += 1
      err = internalnorm(udiff,ttmp)
    end
    zero_vec!(aux_state(utmp))
    state(utmp) .= state(ducache)
    pold_ptr = ptmp

    ldiv!(yprev,Hhalfdt,yprev)
    ldiv!(state(k1),Hhalfdt,state(k1))

    @.. k1 = (k1-utmp)/(dt*ã11)  # r1(y,t) - B1T*z

    ## Stage 2
    k2 = _ode_full_rhs(f,u,pold_ptr,ttmp)
    stats_field(integrator).nf += 1
    @.. k2 *= dt*ã22
    @.. utmp = uprev + k2 + dt*ã21*k1
    ttmp = t + dt*c̃2

    # if applicable, update p, construct new saddle system here, using Hhalfdt
    err, numiter = init_err, init_iter
    u .= utmp
    while err > tol && numiter <= maxiter
      udiff .= u
      ptmp = param_update_func(u,pold_ptr,ttmp)
      S = SaddleSystem(Hhalfdt,f,ptmp,pold_ptr,ducache,solverType;cfact=1.0/(ã22*dt))
      constraint(utmp) .= constraint(_constraint_r2(f,u,ptmp,ttmp))
      ldiv!(mainvector(u),S,mainvector(utmp))
      B1_times_z!(ducache,S)
      @.. udiff -= u
      numiter += 1
      err = internalnorm(udiff,ttmp)
    end
    zero_vec!(aux_state(utmp))
    state(utmp) .= state(ducache)
    pold_ptr = ptmp

    ldiv!(yprev,Hhalfdt,yprev)
    ldiv!(state(k1),Hhalfdt,state(k1))
    ldiv!(state(k2),Hhalfdt,state(k2))

    @.. k2 = (k2-utmp)/(dt*ã22)

    ## Stage 3
    k3 = _ode_full_rhs(f,u,pold_ptr,ttmp)
    stats_field(integrator).nf += 1
    @.. k3 *= dt*ã33
    @.. utmp = uprev + k3 + dt*ã32*k2 + dt*ã31*k1
    ttmp = t + dt

    # if applicable, update p, construct new saddle system here, using Hzero (identity)
    err, numiter = init_err, init_iter
    u .= utmp
    while err > tol && numiter <= maxiter
      udiff .= u
      ptmp = param_update_func(u,pold_ptr,ttmp)
      S = SaddleSystem(Hzero,f,ptmp,pold_ptr,ducache,solverType;cfact=1.0/(ã33*dt))
      constraint(utmp) .= constraint(_constraint_r2(f,u,ptmp,ttmp))
      ldiv!(mainvector(u),S,mainvector(utmp))
      @.. udiff -= u
      numiter += 1
      err = internalnorm(udiff,ttmp)
      #println("error = ",err)
    end
    z = constraint(u)
    @.. z /= (dt*ã33)

    # Final steps
    integrator.p = param_update_func(u,ptmp,t+dt)
    k = f.odef(u, integrator.p, t+dt)
    stats_field(integrator).nf += 1

    integrator.fsallast = k
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.u = u

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
    stats_field(integrator).nf += 1

end

@muladd function perform_step!(integrator,cache::IFHEEulerCache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}
    @unpack t,dt,uprev,u,f,p,opts,alg = integrator
    @unpack internalnorm = opts
    @unpack k1,utmp,udiff,dutmp,fsalfirst,Hdt,S,ptmp,k = cache
    @unpack maxiter, tol = alg
    @unpack param_update_func = f

    init_err = float(1)
    #init_iter = ni ? 1 : maxiter
    init_iter = maxiter  # First-order method does not need iteration

    # aliases to the state and constraint parts
    ytmp, ztmp = state(utmp), constraint(utmp)
    z = constraint(u)
    pold_ptr = p
    pnew_ptr = ptmp

    ttmp = t
    u .= uprev

    _ode_full_rhs!(k1,f,u,pold_ptr,ttmp)
    stats_field(integrator).nf += 1
    @.. k1 *= dt
    @.. utmp = uprev + k1
    ttmp = t + dt

    # if applicable, update p, construct new saddle system here, using Hdt
    err, numiter = init_err, init_iter
    u .= utmp
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew_ptr,u,pold_ptr,ttmp)
      S[1] = SaddleSystem(S[1],Hdt,f,pnew_ptr,pold_ptr,cache)
      _constraint_r2!(utmp,f,u,pnew_ptr,t+dt)
      ldiv!(mainvector(u),S[1],mainvector(utmp))
      @.. udiff -= u
      numiter += 1
      err = internalnorm(udiff,ttmp)
      #println("numiter = ",numiter, ", error = ",err)
    end
    @.. z /= dt

    # Final steps
    param_update_func(p,u,pold_ptr,t)
    f.odef(integrator.fsallast, u, p, t+dt)
    stats_field(integrator).nf += 1

    return nothing
end


function initialize!(integrator,cache::IFHEEulerConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f.odef(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.p = integrator.f.param_update_func(integrator.uprev,integrator.p,integrator.t)
  stats_field(integrator).nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::IFHEEulerConstantCache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}
  @unpack t,dt,uprev,f,p,opts,alg = integrator
  @unpack internalnorm = opts
  @unpack maxiter, tol = alg
  @unpack param_update_func = f

  init_err = float(1)
  #init_iter = ni ? 1 : maxiter
  init_iter = maxiter  # First-order method does not need iteration


  # set up some cache variables
  udiff = deepcopy(uprev)
  ducache = deepcopy(uprev)
  ptmp = deepcopy(p)
  L = _fetch_ode_L(f)
  Hdt = exp(L,-dt,state(uprev))
  pold_ptr = p
  pnew_ptr = ptmp

  k1 = _ode_full_rhs(f,uprev,pold_ptr,t)
  stats_field(integrator).nf += 1
  @.. k1 *= dt
  utmp = @.. uprev + k1

  # if applicable, update p, construct new saddle system here, using Hdt
  err, numiter = init_err, init_iter
  u = deepcopy(utmp)
  while err > tol && numiter <= maxiter
    udiff .= u
    pnew_ptr = param_update_func(u,pold_ptr,t+dt)
    S = SaddleSystem(Hdt,f,pnew_ptr,pold_ptr,ducache,solverType;cfact=1.0/dt)
    constraint(utmp) .= constraint(_constraint_r2(f,u,pnew_ptr,t+dt))
    ldiv!(mainvector(u),S,mainvector(utmp))
    @.. udiff -= u
    numiter += 1
    err = internalnorm(udiff,t+dt)
    #println("error = ",err)
  end
  z = constraint(u)
  @.. z /= dt

  # Final steps
  integrator.p = param_update_func(u,pold_ptr,t)
  k = f.odef(u, integrator.p, t+dt)
  stats_field(integrator).nf += 1

  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end

#### Half-explicit Trapezoidal/Adams-Bashforth 2 (HETrapezoidalAB2) ####

function initialize!(integrator,cache::HETrapezoidalAB2Cache)
    @unpack ki,k,fsalfirst = cache

    integrator.fsalfirst = fsalfirst
    integrator.fsallast = k
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = integrator.fsalfirst
    integrator.k[2] = integrator.fsallast
    integrator.f.odef(integrator.fsalfirst, integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
    #_ode_implicit_rhs!(ki,integrator.f,integrator.uprev,integrator.p,integrator.t)
    integrator.f.param_update_func(integrator.p,integrator.uprev,integrator.p,integrator.t)
    stats_field(integrator).nf += 1

end

@muladd function perform_step!(integrator,cache::HETrapezoidalAB2Cache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}
    @unpack t,dt,uprev,u,f,p,opts,alg = integrator
    @unpack internalnorm = opts
    @unpack ki,ke,utmp,udiff,dutmp,fsalfirst,S,ptmp,k,tab,A = cache
    @unpack α̃1,α̃2,β̃1,β̃2 = tab
    @unpack maxiter, tol = alg
    @unpack param_update_func = f

    init_err = float(1)
    #init_iter = ni ? 1 : maxiter
    init_iter = maxiter
    cnt = integrator.iter

    # aliases to the state and constraint parts
    ytmp, ztmp = state(utmp), constraint(utmp)
    z = constraint(u)
    pold_ptr = p
    pnew_ptr = ptmp

    ttmp = t
    u .= uprev


    if cnt == 1
      # Use Euler step to replace AB2, but still use trapezoidal for implicit part

      @.. utmp = uprev

      # If C is not empty, then find the initial constraint to go along
      # with the initial state
      if _isinvertible(S[1].C)
        _constraint_r2!(utmp,f,utmp,pold_ptr,ttmp)
        constraint_from_state!(mainvector(utmp),S[1],mainvector(utmp))
        @.. ztmp /= (α̃1*dt)
      end

      # Calculate the initial ki (time level 0)
      _ode_implicit_rhs!(ki,f,utmp,pold_ptr,ttmp)

      # Calculate the initial ke (time level 0)
      _ode_r1!(ke,f,utmp,pold_ptr,ttmp)
      stats_field(integrator).nf += 1
      @.. utmp += dt*ke


    else
      @.. utmp = uprev + β̃2*dt*ke

      _ode_r1!(ke,f,uprev,pold_ptr,ttmp)
      stats_field(integrator).nf += 1
      @.. utmp += β̃1*dt*ke  # utmp now corresponds to y* and x*

    end
    ttmp = t + dt

    # if applicable, update p, construct new saddle system here
    err, numiter = init_err, init_iter
    _ode_r1imp!(dutmp,f,utmp,pold_ptr,ttmp)
    @.. u = utmp
    @.. utmp += α̃2*dt*ki + α̃1*dt*dutmp
    while err > tol && numiter <= maxiter
      udiff .= u
      param_update_func(pnew_ptr,u,pold_ptr,ttmp)
      S[1] = SaddleSystem(S[1],A,f,pnew_ptr,pnew_ptr,cache)
      _constraint_r2!(utmp,f,u,pnew_ptr,t+dt)
      ldiv!(mainvector(u),S[1],mainvector(utmp))
      @.. udiff -= u
      numiter += 1
      err = internalnorm(udiff,ttmp)
      #println("numiter = ",numiter, ", error = ",err)
    end

    @.. z /= (α̃1*dt)

    # Final steps
    param_update_func(p,u,pold_ptr,t+dt)
    _ode_implicit_rhs!(ki,f,u,p,t+dt)
    f.odef(integrator.fsallast, u, p, t+dt)
    stats_field(integrator).nf += 1

    return nothing
end


function initialize!(integrator,cache::HETrapezoidalAB2ConstantCache)
  integrator.kshortsize = 2
  integrator.k = typeof(integrator.k)(undef, integrator.kshortsize)
  integrator.fsalfirst = integrator.f.odef(integrator.uprev, integrator.p, integrator.t) # Pre-start fsal
  integrator.p = integrator.f.param_update_func(integrator.uprev,integrator.p,integrator.t)
  stats_field(integrator).nf += 1

  # Avoid undefined entries if k is an array of arrays
  integrator.fsallast = zero(integrator.fsalfirst)
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
end

@muladd function perform_step!(integrator,cache::HETrapezoidalAB2ConstantCache{sc,ni,solverType},repeat_step=false) where {sc,ni,solverType}
  @unpack t,dt,uprev,uprev2,f,p,opts,alg = integrator
  @unpack internalnorm = opts
  @unpack maxiter, tol = alg
  @unpack α̃1,α̃2,β̃1,β̃2 = cache
  @unpack param_update_func = f

  init_err = float(1)
  #init_iter = ni ? 1 : maxiter
  init_iter = maxiter
  cnt = integrator.iter

  # set up some cache variables
  udiff = deepcopy(uprev)
  ducache = deepcopy(uprev)
  ptmp = deepcopy(p)
  L = _fetch_ode_L(f)
  A = implicit_operator(_fetch_ode_L(f),α̃1*dt)
  pold_ptr = p
  pnew_ptr = ptmp

  if cnt == 1
    utmp = uprev
    S = SaddleSystem(A,f,pold_ptr,pold_ptr,ducache,solverType;cfact=1.0/(α̃1*dt))
    if _isinvertible(S.C)
      constraint(utmp) .= constraint(_constraint_r2(f,utmp,pold_ptr,t))
      constraint_from_state!(mainvector(utmp),S,mainvector(utmp))
      ztmp = constraint(utmp)
      @.. ztmp /= (α̃1*dt)
    end

    # Calculate the initial ki (time level 0)
    ki = _ode_implicit_rhs(f,utmp,pold_ptr,t)
    ke = _ode_r1(f,uprev,pold_ptr,t)
    stats_field(integrator).nf += 1
    @.. utmp = utmp + dt*ke

  else
    ke = _ode_r1(f,uprev2,pold_ptr,t-dt) #ke at step n-1
    ki = _ode_implicit_rhs(f,uprev,pold_ptr,t)

    utmp = uprev + β̃2*dt*ke

    ke = _ode_r1(f,uprev,pold_ptr,t)
    stats_field(integrator).nf += 1
    @.. utmp = utmp + β̃1*dt*ke  # utmp now corresponds to y* and x*

  end
  @.. uprev2 = uprev


  ducache .= _ode_r1imp(f,utmp,pold_ptr,t+dt)

  # if applicable, update p, construct new saddle system here, using Hdt
  err, numiter = init_err, init_iter
  u = deepcopy(utmp)
  @.. utmp = utmp + α̃2*dt*ki + α̃1*dt*ducache
  while err > tol && numiter <= maxiter
    udiff .= u
    pnew_ptr = param_update_func(u,pold_ptr,t+dt)
    S = SaddleSystem(A,f,pnew_ptr,pnew_ptr,ducache,solverType;cfact=1.0/(α̃1*dt))
    constraint(utmp) .= constraint(_constraint_r2(f,u,pnew_ptr,t+dt))
    ldiv!(mainvector(u),S,mainvector(utmp))
    @.. udiff -= u
    numiter += 1
    err = internalnorm(udiff,t+dt)
    #println("error = ",err)
  end

  z = constraint(u)
  @.. z /= (α̃1*dt)

  # Final steps
  integrator.p = param_update_func(u,pold_ptr,t)
  k = f.odef(u, integrator.p, t+dt)
  stats_field(integrator).nf += 1

  integrator.fsallast = k
  integrator.k[1] = integrator.fsalfirst
  integrator.k[2] = integrator.fsallast
  integrator.u = u
end
