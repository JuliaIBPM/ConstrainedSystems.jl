import ConstrainedSystems: recursivecopy!, needs_iteration, ArrayPartition

@testset "Iteration test" begin
  Δt = 1e-2
  prob1, _, _ = ConstrainedSystems.basic_constrained_problem()
  integrator = ConstrainedSystems.init(prob1, LiskaIFHERK(),dt=Δt)
  @test needs_iteration(integrator.f,integrator.u,integrator.p,integrator.u) == false

  prob2, _, _ = ConstrainedSystems.cartesian_pendulum_problem()
  integrator = ConstrainedSystems.init(prob2, LiskaIFHERK(),dt=Δt)
  @test needs_iteration(integrator.f,integrator.u,integrator.p,integrator.u) == true

end

y = randn(5)
z = randn(3)
x = randn(2)

@testset "Solution structure" begin

  u = solvector(state=y,constraint=z)
  @test state(u) === y && constraint(u) === z
  @test mainvector(u) === u

  u = solvector(state=y,constraint=z,aux_state=x)
  @test state(u) === y
  @test constraint(u) === z
  @test aux_state(u) === x


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

  r1_c! = ConstrainedSystems._complete_r1(r1!,_func_cache=du)
  r1_c!(du,u,nothing,0.0)

  @test state(du) == fill(1.0,length(y))
  @test constraint(du) == fill(0.0,length(z))
  @test aux_state(du) == fill(2.0,length(x))

  function r2!(dz,p,t)
    fill!(dz,3.0)
  end

  fill!(du,0.0)
  r2_c! = ConstrainedSystems._complete_r2(r2!,_func_cache=du)
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
  B1_c! = ConstrainedSystems._complete_B1(B1!,_func_cache=du)
  B1_c!(du,u,p,0.0)
  @test state(du) == fill(-float(length(z)),length(y))
  @test constraint(du) == fill(0.0,length(z))
  @test aux_state(du) == fill(0.0,length(x))

  fill!(du,0.0)
  fill!(u,1.0)
  p = [B1,B2]
  B2_c! = ConstrainedSystems._complete_B2(B2!,_func_cache=du)
  B2_c!(du,u,p,0.0)
  @test state(du) == fill(0.0,length(y))
  @test constraint(du) == fill(-float(length(y)),length(z))
  @test aux_state(du) == fill(0.0,length(x))

end

@testset "Recursive copy" begin

  struct MyStruct{T}
    a :: T
  end

  n = 1000
  a = rand(n,n)
  z = similar(a)

  s1 = MyStruct(a)
  s2 = MyStruct(z)

  recursivecopy!(s2,s1)
  @test s2.a == s1.a
  @test !(s2.a === s1.a)
  @test @allocated(recursivecopy!(s2,s1)) == 16

  ss1 = MyStruct(s1)
  ss2 = MyStruct(MyStruct(z))
  recursivecopy!(ss2,ss1)
  @test ss2.a.a == ss1.a.a
  @test !(ss2.a.a === ss1.a.a)
  @test @allocated(recursivecopy!(s2,s1)) == 16

  t1 = [s1,s1]
  t2 = [MyStruct(z),MyStruct(z)]
  recursivecopy!(t2,t1)
  @test t2[1].a == t1[1].a
  @test t2[2].a == t1[2].a
  @test !(t2[1].a === t1[1].a)
  @test @allocated(recursivecopy!(s2,s1)) == 16


end
