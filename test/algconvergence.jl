using DiffEqDevTools
import ConstrainedSystems: @unpack

dts = 1 ./ 10 .^(1:4)
testTol = 0.49

struct PendulumParams{P,BT1,BT2}
    params :: P
    B₁ᵀ :: BT1
    B₂ :: BT2
end

θ₀ = 0.0
l = 1.0
g = 1.0
y₀ = Float64[l*cos(θ₀),l*sin(θ₀),0,0]
z₀ = Float64[0.0, 0.0]

u₀ = SaddleVector(y₀,z₀)
du = deepcopy(u₀)

params = [l,g]
p₀ = PendulumParams(params,zeros(Float64,4,2),zeros(Float64,2,4))

function ode_rhs!(dy::Vector{Float64},y::Vector{Float64},p,t)
    dy .= 0.0
    dy[1] = y[3]
    dy[2] = y[4]
    dy[4] = -p.params[2]
    return dy
end

constraint_rhs!(dz::Vector{Float64},p,t) = dz .= [0.0,p.params[1]^2]

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
    B₁ᵀ[3,1] = y[1]; B₁ᵀ[4,1] = y[2]; B₁ᵀ[1,2] = y[1]; B₁ᵀ[2,2] = y[2]
    B₂[1,3] = y[1]; B₂[1,4] = y[2]; B₂[2,1] = y[1]; B₂[2,2] = y[2]
    return q
end




@testset "Convergence test" begin

# Pendulum problem in Cartesian coordinates.

# Get superconverged solution from the basic
# problem expressed in theta
function fex(u,p,t)
    du = similar(u)
    du[1] = u[2]
    du[2] = -p^2*sin(u[1])
    du
  end

u0 = [π/2,0.0]
pex = 1.0  # squared frequency
tspan = (0.0,10.0)
probex = ODEProblem(fex,u0,tspan,pex)
solex = solve(probex, Tsit5(), reltol=1e-16, abstol=1e-16);
xexact(t) = sin(solex(t,idxs=1))
yexact(t) = 1.0 - cos(solex(t,idxs=1))

# Set up the problem
f = ConstrainedODEFunction(ode_rhs!,constraint_rhs!,op_constraint_force!,
                            constraint_op!,
                            _func_cache=deepcopy(du),param_update_func=update_p!)
tspan = (0.0,1.0)
p = deepcopy(p₀)
prob = ODEProblem(f,u₀,tspan,p)


# For now, do this quick and dirty until we figure out how to separate
# state from constraint in sol structure
solutions = [solve(prob,LiskaIFHERK(static_constraints=false);dt=dts[i]) for i=1:length(dts)]

l2err = [sqrt(DiffEqDevTools.recursive_mean(map(x -> float(x).^2,_sol[1,:]-xexact.(_sol.t)))) for _sol in solutions]
error = Dict(:l2 => l2err)

𝒪est = Dict((DiffEqDevTools.calc𝒪estimates(p) for p = pairs(error)))

 @test 𝒪est[:l2][1] ≈ 3 atol=testTol

end
