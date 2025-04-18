using OrdinaryDiffEqTsit5

struct ProblemParams{P,BT1,BT2}
    params :: P
    B₁ᵀ :: BT1
    B₂ :: BT2
end

function basic_unconstrained_problem(;tmax=1.0,iip=true)

  ns = 1
  y0 = 1.0
  α = 1.02
  p = α

  ode_rhs!(dy,y,x,p,t) = dy .= p*y
  ode_rhs(y,x,p,t) = p*y

  y₀ = ones(Float64,ns)
  u₀ = solvector(state=y₀)

  if iip
    f = ConstrainedODEFunction(ode_rhs!,_func_cache=u₀)
  else
    f = ConstrainedODEFunction(ode_rhs)
  end
  tspan = (0.0,1.0)
  prob = ODEProblem(f,u₀,tspan,p)

  xexact(t) = exp(α*t)

  return prob, xexact

end

function basic_unconstrained_if_problem(;tmax=1.0,iip=true)

  ns = 1
  y0 = 1.0
  α = 1.02
  p = α

  ode_rhs!(dy,y,x,p,t) = fill!(dy,0.0)
  ode_rhs(y,x,p,t) = zero(y)

  L = α*I

  y₀ = ones(Float64,ns)
  u₀ = solvector(state=y₀)

  if iip
    f = ConstrainedODEFunction(ode_rhs!,L,_func_cache=u₀)
  else
    f = ConstrainedODEFunction(ode_rhs,L)
  end
  tspan = (0.0,1.0)
  prob = ODEProblem(f,u₀,tspan,p)

  xexact(t) = exp(α*t)

  return prob, xexact

end


function basic_constrained_problem(;tmax=1.0,iip=true)

  U0 = 1.0
  g = 1.0
  α = 0.5
  ω = 5
  y₀ = Float64[0,0,U0,0]
  z₀ = Float64[0]

  params = [U0,g,α,ω];

  p₀ = ProblemParams(params,Array{Float64}(undef,4,1),Array{Float64}(undef,1,4));

  u₀ = solvector(state=y₀,constraint=z₀)
  du = deepcopy(u₀)

  function ode_rhs!(dy::Vector{Float64},y::Vector{Float64},x,p,t)
    dy .= 0.0
    dy[1] = y[3]
    dy[2] = y[4]
    dy[4] = -(y[1]-p.params[1]*t)*p.params[2]
    return dy
  end
  ode_rhs(y::Vector{Float64},x,p,t) = ode_rhs!(deepcopy(y₀),y,x,p,t)

  constraint_rhs!(dz::Vector{Float64},x,p,t) = dz .= Float64[p.params[1]]
  constraint_rhs(x,p,t) = constraint_rhs!(deepcopy(z₀),x,p,t)

  function op_constraint_force!(dy::Vector{Float64},z::Vector{Float64},x,p)
    @unpack B₁ᵀ = p
    dy .= B₁ᵀ*z
  end
  op_constraint_force(z::Vector{Float64},x,p) = op_constraint_force!(deepcopy(y₀),z,x,p)

  function constraint_op!(dz::Vector{Float64},y::Vector{Float64},x,p)
    @unpack B₂ = p
    dz .= B₂*y
  end
  constraint_op(y::Vector{Float64},x,p) = constraint_op!(deepcopy(z₀),y,x,p)

  function update_p!(q,u,p,t)
    y, z = state(u), constraint(u)
    @unpack B₁ᵀ, B₂ = q
    B₁ᵀ .= 0
    B₂ .= 0
    B₁ᵀ[3,1] = 1/(1+q.params[3]*sin(q.params[4]*t))
    B₂[1,3] = 1/(1+q.params[3]*sin(q.params[4]*t))
    return q
  end
  update_p(u,p,t) = update_p!(deepcopy(p),u,p,t)

  if iip
    f = ConstrainedODEFunction(ode_rhs!,constraint_rhs!,op_constraint_force!,
                              constraint_op!,_func_cache=deepcopy(du),
                                            param_update_func=update_p!)
  else
    f = ConstrainedODEFunction(ode_rhs,constraint_rhs,op_constraint_force,
                              constraint_op,param_update_func=update_p)
   end

  tspan = (0.0,tmax)
  p = deepcopy(p₀)
  prob = ODEProblem(f,u₀,tspan,p)

  yexact(t) = g*α*U0/ω*(-0.5*t^2 - cos(ω*t)/ω^2 + 1/ω^2)
  xexact(t) = U0*(t - α*cos(ω*t)/ω+α/ω)

  return prob, xexact, yexact
end


function cartesian_pendulum_problem(;tmax=1.0,iip=true)

  θ₀ = π/2
  l = 1.0
  g = 1.0
  y₀ = Float64[l*sin(θ₀),-l*cos(θ₀),0,0]
  z₀ = Float64[0.0, 0.0]

  u₀ = solvector(state=y₀,constraint=z₀)
  du = deepcopy(u₀)

  params = [l,g]
  p₀ = ProblemParams(params,Array{Float64}(undef,4,2),Array{Float64}(undef,2,4))


  function pendulum_rhs!(dy::Vector{Float64},y::Vector{Float64},x,p,t)
    dy .= 0.0
    dy[1] = y[3]
    dy[2] = y[4]
    dy[4] = -p.params[2]
    return dy
  end

  pendulum_rhs(y::Vector{Float64},x,p,t) = pendulum_rhs!(zero(y),y,x,p,t)


  length_constraint_rhs!(dz::Vector{Float64},x,p,t) = dz .= [0.0,p.params[1]^2]
  length_constraint_rhs(x,p,t) = length_constraint_rhs!(zero(z₀),x,p,t)


  function length_constraint_force!(dy::Vector{Float64},z::Vector{Float64},x,p)
    @unpack B₁ᵀ = p
    dy .= B₁ᵀ*z
  end
  length_constraint_force(z::Vector{Float64},x,p) = length_constraint_force!(zero(y₀),z,x,p)

  function length_constraint_op!(dz::Vector{Float64},y::Vector{Float64},x,p)
    @unpack B₂ = p
    dz .= B₂*y
  end
  length_constraint_op(y::Vector{Float64},x,p) = length_constraint_op!(zero(z₀),y,x,p)

  function update_p!(q,u,p,t)
    y, z = state(u), constraint(u)
    @unpack B₁ᵀ, B₂ = q
    fill!(B₁ᵀ,0.0)
    fill!(B₂,0.0)
    B₁ᵀ[3,1] = y[1]; B₁ᵀ[4,1] = y[2]; B₁ᵀ[1,2] = y[1]; B₁ᵀ[2,2] = y[2]
    B₂[1,3] = y[1]; B₂[1,4] = y[2]; B₂[2,1] = y[1]; B₂[2,2] = y[2]
    return q
  end
  update_p(u,p,t) = update_p!(deepcopy(p),u,p,t)

  if iip
    f = ConstrainedODEFunction(pendulum_rhs!,length_constraint_rhs!,length_constraint_force!,
                                length_constraint_op!,
                               _func_cache=deepcopy(du),param_update_func=update_p!)
  else
    f = ConstrainedODEFunction(pendulum_rhs,length_constraint_rhs,length_constraint_force,
                                length_constraint_op,param_update_func=update_p)
  end

  tspan = (0.0,tmax)
  p = deepcopy(p₀)
  prob = ODEProblem(f,u₀,tspan,p)

  # Get superconverged solution from the basic
  # problem expressed in theta
  function pendulum_theta(u,p,t)
      du = similar(u)
      du[1] = u[2]
      du[2] = -p^2*sin(u[1])
      du
  end

  u0 = [π/2,0.0]
  pex = g/l  # squared frequency
  tspan = (0.0,10.0)
  probex = ODEProblem(pendulum_theta,u0,tspan,pex)
  solex = solve(probex, Tsit5(), reltol=1e-16, abstol=1e-16);
  xexact(t) = sin(solex(t,idxs=1))
  yexact(t) = -cos(solex(t,idxs=1))

  return prob, xexact, yexact

end

function partitioned_problem(;tmax=1.0,iip=true)

  ω = 1.0
  βu = -0.2
  βv = -0.5

  par = [ω,βu,βv]

  U₀ = Float64[0,1]
  X₀ = Float64[1,0]
  Z₀ = Float64[0]
  u₀ = solvector(state=X₀,constraint=Z₀,aux_state=U₀)
  du = deepcopy(u₀)

  B₂ = Array{Float64}(undef,1,2)
  B₁ᵀ = Array{Float64}(undef,2,1)

  p₀ = ProblemParams(par,B₁ᵀ,B₂)

  L = Diagonal([βu,βv])

  function X_rhs!(dy,y,x,p,t)
    fill!(dy,0.0)
    return dy
  end
  X_rhs(y,x,p,t) = X_rhs!(deepcopy(X₀),y,x,p,t)

  function U_rhs!(dy,u,p,t)
    fill!(dy,0.0)
    ω = p.params[1]
    dy[1] = -ω*cos(ω*t)
    dy[2] = -ω*sin(ω*t)
    return dy
  end
  U_rhs(u,p,t) = U_rhs!(deepcopy(U₀),u,p,t)

  ode_rhs! = ArrayPartition((X_rhs!,U_rhs!))
  ode_rhs = ArrayPartition((X_rhs,U_rhs))

  constraint_rhs!(dz,x,p,t) = dz .= Float64[0]
  constraint_rhs(x,p,t) = constraint_rhs!(deepcopy(Z₀),x,p,t)

  function op_constraint_force!(dy,z,x,p)
    @unpack B₁ᵀ = p
    dy .= B₁ᵀ*z
    return dy
  end
  op_constraint_force(z,x,p) = op_constraint_force!(deepcopy(X₀),z,x,p)

  function constraint_op!(dz,y,x,p)
    @unpack B₂ = p
    dz .= B₂*y
  end
  constraint_op(y,x,p) = constraint_op!(deepcopy(Z₀),y,x,p)

  function update_p!(q,u,p,t)
    x = aux_state(u)
    @unpack B₁ᵀ, B₂ = q

    B₁ᵀ[1,1] = x[1]
    B₁ᵀ[2,1] = x[2]
    B₂[1,1] = x[1]
    B₂[1,2] = x[2]
    return q
  end
  update_p(u,p,t) = update_p!(deepcopy(p),u,p,t)

  if iip
    f = ConstrainedODEFunction(ode_rhs!,constraint_rhs!,op_constraint_force!,
                              constraint_op!,L,_func_cache=deepcopy(du),
                              param_update_func=update_p!)
  else
    f = ConstrainedODEFunction(ode_rhs,constraint_rhs,op_constraint_force,
                              constraint_op,L,param_update_func=update_p)
  end

  tspan = (0.0,tmax)
  p = deepcopy(p₀)
  update_p!(p,u₀,p,0.0)
  prob = ODEProblem(f,u₀,tspan,p)

  # function fex(du,u,p,t)
  #   UV = u.x[1]
  #   xy = u.x[2]
  #   ω = p[1]
  #   βu = p[2]
  #   βv = p[3]
  #   du[1] = -ω*cos(ω*t)
  #   du[2] = -ω*sin(ω*t)
  #   Usq = u[1]^2+u[2]^2
  #   du[3] = (βu-u[1]/Usq*(du[1]+βu*u[1]))*u[3] - u[1]/Usq*(du[2]+βv*u[2])*u[4]
  #   du[4] = -u[2]/Usq*(du[1]+βu*u[1])*u[3] + (βv-u[2]/Usq*(du[2]+βv*u[2]))*u[4]
  #   return nothing
  # end
  #
  # tspan = (0.0,tmax)
  # u₀ex = SaddleVector(U₀,X₀)
  #
  # probex = ODEProblem(fex,u₀ex,tspan,par)
  # solex = solve(probex, Tsit5(), reltol=1e-16, abstol=1e-16)
  # xexact(t) = solex(t,idxs=3)
  # yexact(t) = solex(t,idxs=4)


  fex(t) = exp(0.5*(βu+βv)*t)*exp(0.25/ω*sin(2ω*t)*(βu-βv))
  xexact(t) = cos(ω*t)*fex(t)
  yexact(t) = sin(ω*t)*fex(t)


  return prob, xexact, yexact

end

function basic_constrained_if_problem_with_cmatrix(;tmax=1.0,iip=true)

  y0 = 1
  y₀ = Float64[y0]
  z₀ = Float64[0]
  u₀ = solvector(state=y₀,constraint=z₀)
  du = deepcopy(u₀)

  α = -1.0
  B1T = 1.0
  B2 = 1.0
  β = 1.2
  r2 = 1.0
  p = [α,B1T,B2,β,r2]

  L = α*I(1)

  ode_rhs!(dy,y,x,p,t) = fill!(dy,0.0)
  ode_rhs(y,x,p,t) = zero(y)

  constraint_rhs!(dz,x,p,t) = fill!(dz,p[5])
  constraint_rhs(x,p,t) = p[5]*ones(1)

  function op_constraint_force!(dy::Vector{Float64},z::Vector{Float64},x,p)
    dy .= p[2]*z
  end

  op_constraint_force(z::Vector{Float64},x,p) = op_constraint_force!(deepcopy(y₀),z,x,p)

  function constraint_op!(dz::Vector{Float64},y::Vector{Float64},x,p)
    dz .= p[3]*y
  end
  constraint_op(y::Vector{Float64},x,p) = constraint_op!(deepcopy(z₀),y,x,p)

  function constraint_reg!(dz::Vector{Float64},z::Vector{Float64},x,p)
    dz .= p[4]*z
  end
  constraint_reg(z::Vector{Float64},x,p) = constraint_reg!(deepcopy(z₀),z,x,p)

  if iip
    f = ConstrainedODEFunction(ode_rhs!,constraint_rhs!,op_constraint_force!,
                              constraint_op!,L,constraint_reg!,_func_cache=deepcopy(du))
  else
    f = ConstrainedODEFunction(ode_rhs,constraint_rhs,op_constraint_force,
                              constraint_op,L,constraint_reg)
   end

  tspan = (0.0,tmax)
  prob = ODEProblem(f,u₀,tspan,p)

  xexact(t) = exp((α+1/β)*t)*y0 + r2/β/(α+1/β)*(1 - exp((α+1/β)*t))
  yexact(t) = (r2 - xexact(t))/β

  return prob, xexact, yexact
end
