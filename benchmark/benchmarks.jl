using BenchmarkTools, ConstrainedSystems

Ns, Nc = 1000, 100
A,B2,B1t,C,rhs1v,rhs2v,solex = ConstrainedSystems.basic_linalg_problem(Ns=Ns,Nc=Nc);

Δt = 1.0e-2
prob1, xexact, yexact = ConstrainedSystems.basic_constrained_problem(tmax=3Δt)
prob2, xexact, yexact = ConstrainedSystems.cartesian_pendulum_problem(tmax=3Δt)


BenchmarkTools.DEFAULT_PARAMETERS.gcsample = true

SUITE = BenchmarkGroup()

SUITE["saddle setup"] = BenchmarkGroup()
SUITE["saddle solve"] = BenchmarkGroup()

rhs = SaddleVector(rhs1v,rhs2v)
sol = similar(rhs)

SUITE["saddle setup"]["basic problem"] = @benchmarkable ConstrainedSystems.SaddleSystem(A,B2,B1t,C,rhs)

As = ConstrainedSystems.SaddleSystem(A,B2,B1t,C,rhs)

SUITE["saddle solve"]["basic problem"] = @benchmarkable sol .= As\rhs


SUITE["timemarch setup"] = BenchmarkGroup()
SUITE["timemarch basic problem"] = BenchmarkGroup()
SUITE["timemarch pendulum problem"] = BenchmarkGroup()

SUITE["timemarch setup"]["IFHEEuler"] = @benchmarkable ConstrainedSystems.init(prob2, IFHEEuler(),dt=Δt)
SUITE["timemarch setup"]["LiskaIFHERK"] = @benchmarkable ConstrainedSystems.init(prob2, LiskaIFHERK(),dt=Δt)

integrator11 = ConstrainedSystems.init(prob1, IFHEEuler(),dt=Δt)
integrator12 = ConstrainedSystems.init(prob1, LiskaIFHERK(),dt=Δt)
integrator21 = ConstrainedSystems.init(prob2, IFHEEuler(),dt=Δt)
integrator22 = ConstrainedSystems.init(prob2, LiskaIFHERK(),dt=Δt)

#SUITE["basic problem"]["LiskaIFHERK"] = @benchmarkable solve(prob1, LiskaIFHERK(),dt=Δt)
#SUITE["basic problem"]["IFHEEuler"] = @benchmarkable solve(prob1, IFHEEuler(),dt=Δt)
#SUITE["pendulum problem"]["LiskaIFHERK"] = @benchmarkable solve(prob2, LiskaIFHERK(),dt=Δt)
#SUITE["pendulum problem"]["IFHEEuler"] = @benchmarkable solve(prob2, IFHEEuler(),dt=Δt)

SUITE["timemarch basic problem"]["IFHEEuler"] = @benchmarkable step!(integrator11,Δt)
SUITE["timemarch basic problem"]["LiskaIFHERK"] = @benchmarkable step!(integrator12,Δt)
SUITE["timemarch pendulum problem"]["IFHEEuler"] = @benchmarkable step!(integrator21,Δt)
SUITE["timemarch pendulum problem"]["LiskaIFHERK"] = @benchmarkable step!(integrator22,Δt)

run(SUITE)
