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

end
