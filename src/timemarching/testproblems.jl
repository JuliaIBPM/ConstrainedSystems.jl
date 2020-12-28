struct ProblemParams{P,BT1,BT2}
    params :: P
    B₁ᵀ :: BT1
    B₂ :: BT2
end


function basic_constrained_problem(;tmax=1.0)

  U0 = 1.0
  g = 1.0
  α = 0.5
  ω = 5
  y₀ = Float64[0,0,U0,0]
  z₀ = Float64[0]

  params = [U0,g,α,ω];

  p₀ = ProblemParams(params,Array{Float64}(undef,4,1),Array{Float64}(undef,1,4));

  u₀ = SaddleVector(y₀,z₀)
  du = deepcopy(u₀)

  function ode_rhs!(dy::Vector{Float64},y::Vector{Float64},p,t)
    dy .= 0.0
    dy[1] = y[3]
    dy[2] = y[4]
    dy[4] = -(y[1]-p.params[1]*t)*p.params[2]
    return dy
  end

  constraint_rhs!(dz::Vector{Float64},p,t) = dz .= Float64[p.params[1]]

  function op_constraint_force!(dy::Vector{Float64},z::Vector{Float64},p)
    @unpack B₁ᵀ = p
    dy .= B₁ᵀ*z
  end

  function constraint_op!(dz::Vector{Float64},y::Vector{Float64},p)
    @unpack B₂ = p
    dz .= B₂*y
  end

  function update_p!(q,u,p,t)
    y, z = state(u), constraint(u)
    @unpack B₁ᵀ, B₂ = q
    B₁ᵀ .= 0
    B₂ .= 0
    B₁ᵀ[3,1] = 1/(1+q.params[3]*sin(q.params[4]*t))
    B₂[1,3] = 1/(1+q.params[3]*sin(q.params[4]*t))
    return q
  end

  f = ConstrainedODEFunction(ode_rhs!,constraint_rhs!,op_constraint_force!,
                              constraint_op!,_func_cache=deepcopy(du),
                                            param_update_func=update_p!)

  tspan = (0.0,tmax)
  p = deepcopy(p₀)
  prob = ODEProblem(f,u₀,tspan,p)

  yexact(t) = g*α*U0/ω*(-0.5*t^2 - cos(ω*t)/ω^2 + 1/ω^2)
  xexact(t) = U0*(t - α*cos(ω*t)/ω+α/ω)

  return prob, xexact, yexact
end


function cartesian_pendulum_problem(;tmax=1.0)

  θ₀ = π/2
  l = 1.0
  g = 1.0
  y₀ = Float64[l*sin(θ₀),-l*cos(θ₀),0,0]
  z₀ = Float64[0.0, 0.0]

  u₀ = SaddleVector(y₀,z₀)
  du = deepcopy(u₀)

  params = [l,g]
  p₀ = ProblemParams(params,Array{Float64}(undef,4,2),Array{Float64}(undef,2,4))

  function pendulum_rhs!(dy::Vector{Float64},y::Vector{Float64},p,t)
    dy .= 0.0
    dy[1] = y[3]
    dy[2] = y[4]
    dy[4] = -p.params[2]
    return dy
  end

  length_constraint_rhs!(dz::Vector{Float64},p,t) = dz .= [0.0,p.params[1]^2]

  function length_constraint_force!(dy::Vector{Float64},z::Vector{Float64},p)
    @unpack B₁ᵀ = p
    dy .= B₁ᵀ*z
  end

  function length_constraint_op!(dz::Vector{Float64},y::Vector{Float64},p)
    @unpack B₂ = p
    dz .= B₂*y
  end

  function update_p!(q,u,p,t)
    y, z = state(u), constraint(u)
    @unpack B₁ᵀ, B₂ = q
    fill!(B₁ᵀ,0.0)
    fill!(B₂,0.0)
    B₁ᵀ[3,1] = y[1]; B₁ᵀ[4,1] = y[2]; B₁ᵀ[1,2] = y[1]; B₁ᵀ[2,2] = y[2]
    B₂[1,3] = y[1]; B₂[1,4] = y[2]; B₂[2,1] = y[1]; B₂[2,2] = y[2]
    return q
  end

  f = ConstrainedODEFunction(pendulum_rhs!,length_constraint_rhs!,length_constraint_force!,
                              length_constraint_op!,
                              _func_cache=deepcopy(du),param_update_func=update_p!)
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

function partitioned_problem(;tmax=1.0)

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

  function X_rhs!(dy,y,p,t)
    fill!(dy,0.0)
    return dy
  end

  function U_rhs!(dy,y,p,t)
    fill!(dy,0.0)
    ω = p.params[1]
    dy[1] = -ω*cos(ω*t)
    dy[2] = -ω*sin(ω*t)
    return dy
  end

  ode_rhs! = ArrayPartition((X_rhs!,U_rhs!))

  constraint_rhs!(dz,p,t) = dz .= Float64[0]


  function op_constraint_force!(dy,z,p)
    @unpack B₁ᵀ = p
    dy .= B₁ᵀ*z
    return dy
  end

  function constraint_op!(dz,y,p)
    @unpack B₂ = p
    dz .= B₂*y
  end

  function update_p!(q,u,p,t)
    x = aux_state(u)
    @unpack B₁ᵀ, B₂ = q

    B₁ᵀ[1,1] = x[1]
    B₁ᵀ[2,1] = x[2]
    B₂[1,1] = x[1]
    B₂[1,2] = x[2]
    return q
  end

  f = ConstrainedODEFunction(ode_rhs!,constraint_rhs!,op_constraint_force!,
                              constraint_op!,L,_func_cache=deepcopy(du),
                              param_update_func=update_p!)
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

# out of place

function basic_constrained_problem_oop(;tmax=1.0)

  U0 = 1.0
  g = 1.0
  α = 0.5
  ω = 5
  y₀ = Float64[0,0,U0,0]
  z₀ = Float64[0]

  params = [U0,g,α,ω];

  p₀ = ProblemParams(params,Array{Float64}(undef,4,1),Array{Float64}(undef,1,4));

  u₀ = SaddleVector(y₀,z₀)
  du = deepcopy(u₀)

  function ode_rhs(y::Vector{Float64},p,t)
    dy = zero(y)
    dy[1] = y[3]
    dy[2] = y[4]
    dy[4] = -(y[1]-p.params[1]*t)*p.params[2]
    return dy
  end

  constraint_rhs(p,t) = Float64[p.params[1]]

  op_constraint_force(z::Vector{Float64},p) = p.B₁ᵀ*z

  constraint_op(y::Vector{Float64},p) = p.B₂*y

  function update_p(u,p,t)
    q = deepcopy(p)

    @unpack B₁ᵀ, B₂ = q
    B₁ᵀ .= 0
    B₂ .= 0
    B₁ᵀ[3,1] = 1/(1+q.params[3]*sin(q.params[4]*t))
    B₂[1,3] = 1/(1+q.params[3]*sin(q.params[4]*t))
    return q
  end

  f = ConstrainedODEFunction(ode_rhs,constraint_rhs,op_constraint_force,
                              constraint_op,param_update_func=update_p)

  tspan = (0.0,tmax)
  p = deepcopy(p₀)
  prob = ODEProblem(f,u₀,tspan,p)

  yexact(t) = g*α*U0/ω*(-0.5*t^2 - cos(ω*t)/ω^2 + 1/ω^2)
  xexact(t) = U0*(t - α*cos(ω*t)/ω+α/ω)

  return prob, xexact, yexact
end
