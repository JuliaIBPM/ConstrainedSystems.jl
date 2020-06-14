using CartesianGrids
using LinearAlgebra

@testset "Time Marching" begin

  @testset "Basic Unconstrained" begin

  ω = 4
  u₀ = 1.0
  uex(t) = u₀ + sin(ω*t)/ω

  Δt = 0.005
  T = 0:Δt:10
  u = [u₀]
  r₁(u::Vector{Float64},t::Float64) = cos(ω*t)
  rk = RK(u,Δt,r₁,rk=ConstrainedSystems.RK31)

  u = [u₀]
  uhist = Float64[]
  for t in T
    push!(uhist,u[1])
    t,u = rk(t,u)
  end

  @test norm(uhist-uex.(T)) ≈ 0.000494857145025

  end

  @testset "Unconstrained with integrating factor" begin

  α = 0.5
  ω = 4
  u₀ = 1.0
  uex(t) = u₀*exp(-α*t) + (α*(cos(ω*t)-exp(-α*t))+ω*sin(ω*t))/(α^2+ω^2)

  CartesianGrids.plan_intfact(t::Float64,u::Vector{Float64}) = exp(-α*t)

  Δt = 0.005
  T = 0:Δt:10
  u = [u₀]
  r₁(u::Vector{Float64},t::Float64) = cos(ω*t)
  ifrk = IFRK(u,Δt,plan_intfact,r₁,rk=ConstrainedSystems.RK31)

  u = [u₀]
  uhist = Float64[]
  for t in T
    push!(uhist,u[1])
    t,u = ifrk(t,u)
  end
  @test norm(uhist-uex.(T)) ≈ 0.0005014844449

  end

  @testset "Basic unconstrained using IFHERK" begin

  ω = 4
  u₀ = 1.0
  uex(t) = u₀ + sin(ω*t)/ω

  Δt = 0.005
  T = 0:Δt:10
  u = [u₀]
  f = Vector{Float64}()
  r₁(u::Vector{Float64},t::Float64) = cos(ω*t)
  r₂(u::Vector{Float64},t::Float64) = Vector{Float64}()
  plan_constraints(u::Vector{Float64},t::Float64) = f -> zeros(Float64,1), u -> Vector{Float64}()
  CartesianGrids.plan_intfact(t::Float64,u::Vector{Float64}) = Matrix(1.0I,1,1)
  ifherk = IFHERK(u,f,Δt,plan_intfact,plan_constraints,(r₁,r₂),rk=ConstrainedSystems.RK31)

  u = [u₀]
  uhist = Float64[]
  for t in T
    push!(uhist,u[1])
    t,u,_ = ifherk(t,u)
  end

  @test norm(uhist-uex.(T)) ≈ 0.0004948571450253

  end

  @testset "Constrained using IFHERK" begin

  nx = 129; ny = 129; Lx = 2.0; Δx = Lx/(nx-2)
  u₀ = Nodes(Dual,(nx,ny)); # field initial condition

  n = 128; θ = range(0,stop=2π,length=n+1)
  R = 0.5; xb = 1.0 .+ R*cos.(θ); yb = 1.0 .+ R*sin.(θ)
  X = VectorData(xb[1:n],yb[1:n])
  f = ScalarData(X); # to be used as the Lagrange multiplier

  reg = Regularize(X,Δx;issymmetric=true)
  Hmat, Emat = RegularizationMatrix(reg,f,u₀)
  plan_constraints(u::Nodes{Dual,nx,ny},t::Float64) = Hmat, Emat

  r₁(u::Nodes{T,NX,NY},t::Float64) where {T,NX,NY} = Nodes(T,u); # sets to zeros
  r₂(u::Nodes{T,NX,NY},t::Float64) where {T,NX,NY} = 1.0; # sets uniformly to 1.0

  Δt = 1.0
  t = 0.0
  u = deepcopy(u₀)

  solver = IFHERK(u,f,Δt,CartesianGrids.plan_intfact,plan_constraints,(r₁,r₂),
                rk=ConstrainedSystems.RK31)


  for i = 1:20
      t, u, f = solver(t,u)
  end

  @test isapprox(LinearAlgebra.norm(u),0.3008,atol=1e-4)

  end

end
