using DiffEqDevTools
import DiffEqDevTools: recursive_mean
import ConstrainedSystems: @unpack

dts = 1 ./ 2 .^(9:-1:5)
testTol = 0.2

@inline function compute_l2err(sol,t,sol_analytic)
  sqrt(recursive_mean(map(x -> float(x).^2,sol-sol_analytic.(t))))
end

function compute𝒪est(solutions,idx,sol_analytic)
  l2err = [compute_l2err(_sol[idx,:],_sol.t,sol_analytic) for _sol in solutions]
  error = Dict(:l2 => l2err)
  𝒪est = Dict((DiffEqDevTools.calc𝒪estimates(p) for p = pairs(error)))
end


@testset "Convergence test" begin

prob, xexact, yexact = ConstrainedSystems.basic_constrained_problem()

# For now, do this quick and dirty until we figure out how to separate
# state from constraint in sol structure
solutions1 = [solve(prob,LiskaIFHERK();dt=dts[i]) for i=1:length(dts)]
solutions2 = [solve(prob,IFHEEuler();dt=dts[i]) for i=1:length(dts)]

𝒪est1 = compute𝒪est(solutions1,1,xexact)
𝒪est2 = compute𝒪est(solutions2,1,xexact)

@test 𝒪est1[:l2][1] ≈ 2 atol=testTol
@test 𝒪est2[:l2][1] ≈ 1 atol=testTol


end
