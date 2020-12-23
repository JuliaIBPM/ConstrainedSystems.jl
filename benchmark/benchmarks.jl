using BenchmarkTools, ConstrainedSystems

Δt = 1.0e-2
prob1, xexact, yexact = ConstrainedSystems.basic_constrained_problem(tmax=3Δt)
prob2, xexact, yexact = ConstrainedSystems.cartesian_pendulum_problem(tmax=3Δt)

BenchmarkTools.DEFAULT_PARAMETERS.gcsample = true

SUITE = BenchmarkGroup()
SUITE["basic problem"] = BenchmarkGroup()
SUITE["pendulum problem"] = BenchmarkGroup()

SUITE["basic problem"]["LiskaIFHERK"] = @benchmarkable solve(prob1, LiskaIFHERK(),dt=Δt)
SUITE["basic problem"]["IFHEEuler"] = @benchmarkable solve(prob1, IFHEEuler(),dt=Δt)
SUITE["pendulum problem"]["LiskaIFHERK"] = @benchmarkable solve(prob2, LiskaIFHERK(),dt=Δt)
SUITE["pendulum problem"]["IFHEEuler"] = @benchmarkable solve(prob2, IFHEEuler(),dt=Δt)

run(SUITE)
