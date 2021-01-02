using LinearAlgebra

ns = 5
nc = 3
na = 2
y = randn(ns)
z = randn(nc)
x = randn(na)

@testset "Solution structure" begin

  u = solvector(state=y,constraint=z)
  @test state(u) === y && constraint(u) === z
  @test mainvector(u) === ArrayPartition(y,z)

  u = solvector(state=y,constraint=z,aux_state=x)
  @test state(u) === y
  @test constraint(u) === z
  @test aux_state(u) === x

  u = solvector(state=y)
  @test state(u) === y
  @test constraint(u) == empty(y)
  @test aux_state(u) == nothing
  @test mainvector(u) == ArrayPartition(y,empty(y))

  e = _empty(y)
  @test isempty(e)

end

@testset "Function structure" begin

  u = solvector(state=y,constraint=z,aux_state=x)
  du = similar(u)
  fill!(du,0.0)

  function state_r1!(dy,y,p,t)
    fill!(dy,1.0)
  end

  function aux_r1!(dy,y,p,t)
    fill!(dy,2.0)
  end

  r1! = ArrayPartition((state_r1!,aux_r1!))

  r1_c! = ConstrainedSystems._complete_r1(r1!,Val(true),_func_cache=du)
  r1_c!(du,u,nothing,0.0)

  @test state(du) == fill(1.0,length(y))
  @test constraint(du) == fill(0.0,length(z))
  @test aux_state(du) == fill(2.0,length(x))

  function r2!(dz,p,t)
    fill!(dz,3.0)
  end

  fill!(du,0.0)
  r2_c! = ConstrainedSystems._complete_r2(r2!,Val(true),_func_cache=du)
  r2_c!(du,u,nothing,0.0)
  @test state(du) == fill(0.0,length(y))
  @test constraint(du) == fill(3.0,length(z))
  @test aux_state(du) == fill(0.0,length(x))

  B1 = ones(Float64,length(y),length(z))
  B2 = ones(Float64,length(z),length(y))

  function B1!(dy,z,p)
    dy .= p[1]*z
  end
  function B2!(dz,y,p)
    dz .= p[2]*y
  end

  fill!(du,0.0)
  fill!(u,1.0)
  p = [B1,B2]
  B1_c! = ConstrainedSystems._complete_B1(B1!,Val(true),_func_cache=du)
  B1_c!(du,u,p,0.0)
  @test state(du) == fill(-float(length(z)),length(y))
  @test constraint(du) == fill(0.0,length(z))
  @test aux_state(du) == fill(0.0,length(x))

  fill!(du,0.0)
  fill!(u,1.0)
  p = [B1,B2]
  B2_c! = ConstrainedSystems._complete_B2(B2!,Val(true),_func_cache=du)
  B2_c!(du,u,p,0.0)
  @test state(du) == fill(0.0,length(y))
  @test constraint(du) == fill(-float(length(y)),length(z))
  @test aux_state(du) == fill(0.0,length(x))

end

@testset "ConstrainedODEFunction" begin

  ns = 5
  nc = 2
  na = 3
  y = ones(Float64,ns)
  z = ones(Float64,nc)
  x = ones(Float64,na)

  B1 = Array{Float64}(undef,ns,nc)
  B2 = Array{Float64}(undef,nc,ns)
  B1 .= 1:ns
  B2 .= transpose(B1)
  p = [B1,B2]

  ode_rhs!(dy,y,p,t) = dy .= 1.01*y
  constraint_force!(dy,z,p) = dy .= p[1]*z
  constraint_rhs!(dz,p,t) = dz .= 1.0
  constraint_op!(dz,y,p) = dz .= p[2]*y

  u₀ = solvector(state=y,constraint=z)
  du = zero(u₀)

  f = ConstrainedODEFunction(ode_rhs!,constraint_rhs!,
                              constraint_force!,constraint_op!,_func_cache=u₀)
  f(du,u₀,p,0.0)
  @test state(du) == 1.01*y .- B1*z
  @test constraint(du) == 1.0 .- B2*y

  ode_rhs(y,p,t) = 1.01*y
  constraint_force(z,p) = p[1]*z
  constraint_rhs(p,t) = fill(1.0,nc)
  constraint_op(y,p) = p[2]*y

  f = ConstrainedODEFunction(ode_rhs,constraint_rhs,
                            constraint_force,constraint_op)
  du = f(u₀,p,0.0)
  @test state(du) == 1.01*y .- B1*z
  @test constraint(du) == 1.0 .- B2*y

  L = 2*I
  f = ConstrainedODEFunction(ode_rhs!,constraint_rhs!,
                              constraint_force!,constraint_op!,L,_func_cache=u₀)

  f(du,u₀,p,0.0)
  @test state(du) == 1.01*y .- B1*z .+ L*y
  @test state(du) ≈ [1.01,-0.99,-2.99,-4.99,-6.99] atol=1e-12
  @test constraint(du) == 1.0 .- B2*y
  @test constraint(du) ≈ [-14.0, -14.0] atol=1e-12

  f = ConstrainedODEFunction(ode_rhs,constraint_rhs,
                            constraint_force,constraint_op,L)

  du = f(u₀,p,0.0)
  @test state(du) == 1.01*y .- B1*z .+ L*y
  @test state(du) ≈ [1.01,-0.99,-2.99,-4.99,-6.99] atol=1e-12
  @test constraint(du) == 1.0 .- B2*y
  @test constraint(du) ≈ [-14.0, -14.0] atol=1e-12

  # Unconstrained system
  u₀ = solvector(state=y)
  du = zero(u₀)
  f = ConstrainedODEFunction(ode_rhs!,L,_func_cache=u₀)
  f(du,u₀,p,0.0)
  @test state(du) == 1.01*y .+ L*y
  @test state(du) ≈ fill(3.01,ns) atol=1e-12
  @test constraint(du) == empty(y)

  # Test that unconstrained systems create a saddle system that works properly
  # Should only invert the upper left operator
  u = deepcopy(u₀)
  S = SaddleSystem(2*I,f,p,p,deepcopy(du),Direct)
  u .= S\u₀
  @test state(u) ≈ fill(0.5,ns) atol=1e-12
  @test constraint(u) == empty(y)

  f = ConstrainedODEFunction(ode_rhs,L)
  du = f(u₀,p,0.0)
  @test state(du) == 1.01*y .+ L*y
  @test state(du) ≈ fill(3.01,ns) atol=1e-12
  @test constraint(du) == empty(y)

  # Test that unconstrained systems create a saddle system that works properly
  # Should only invert the upper left operator
  u = deepcopy(u₀)
  S = SaddleSystem(2*I,f,p,p,deepcopy(du),Direct)
  u .= S\u₀
  @test state(u) ≈ fill(0.5,ns) atol=1e-12
  @test constraint(u) == empty(y)


end
