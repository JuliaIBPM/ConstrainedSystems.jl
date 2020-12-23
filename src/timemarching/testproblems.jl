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
