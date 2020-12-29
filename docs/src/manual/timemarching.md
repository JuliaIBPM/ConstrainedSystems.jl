# Time marching

```@meta
DocTestSetup = quote
using ConstrainedSystems
using CartesianGrids
end
```

```math
\def\ddt#1{\frac{\mathrm{d}#1}{\mathrm{d}t}}

\renewcommand{\vec}{\boldsymbol}
\newcommand{\uvec}[1]{\vec{\hat{#1}}}
\newcommand{\utangent}{\uvec{\tau}}
\newcommand{\unormal}{\uvec{n}}

\renewcommand{\d}{\,\mathrm{d}}
```


```@setup create
using ConstrainedSystems
using CartesianGrids
using Plots
```
`ConstrainedSystems` is equipped with tools for solving systems of equations of the
general form of half-explicit differential-algebraic equations,

$$\ddt y = L u - B_1^T(y,t) z + r_1(y,t), \quad B_2(y,t) y = r_2(t), \quad y(0) = y_0$$

where $z$ is the Lagrange multiplier for enforcing the constraints on $y$. Note
that the constraint operators may depend on the state and on time. The linear operator $L$ may be a matrix or a scalar, but is generally independent of time. (The method of integrating factors can deal with time-dependent $L$, but we don't encounter such systems in the `ConstrainedSystems` context so we won't discuss them.) Our objective is to solve
for $y(t)$ and $z(t)$.


## Constrained integrating factor systems


Let's demonstrate this on the example of heat diffusion from a circular ring whose temperature
is held constant. In this case, $L$ is the discrete Laplace operator times the heat diffusivity,
$r_1$ is zero (in the absence of volumetric heating sources), and $r_2$ is the temperature of
the ring. The operators $B_1^T$ and $B_2$ will be the regularization and interpolation
operators between discrete point-wise data on the ring and the field data.

The ring will have radius $1/2$ and fixed temperature $1$, and
the heat diffusivity is $1$. (In other words, the problem has been non-dimensionalized
by the diameter of the circle, the dimensional ring temperature, and the dimensional diffusivity.)

First, we will construct a field to accept the temperature on

```@repl march
nx = 129; ny = 129; Lx = 2.0; Δx = Lx/(nx-2);
w₀ = Nodes(Dual,(nx,ny)); # field initial condition
```

Now set up a ring of points on the circle at center $(1,1)$.

```@repl march
n = 128; θ = range(0,stop=2π,length=n+1);
R = 0.5; xb = 1.0 .+ R*cos.(θ); yb = 1.0 .+ R*sin.(θ);
X = VectorData(xb[1:n],yb[1:n]);
z = ScalarData(X); # to be used as the Lagrange multiplier
```

Together, `w₀` and `z` comprise the initial solution vector:

```@repl march
u₀ = solvector(state=w₀,constraint=z)
```

Now set up the operators. We first set up the linear operator, a Laplacian endowed
with its inverse:

```@repl march
L = plan_laplacian(w₀,with_inverse=true)
```

Now the right-hand side operators for the ODEs and constraints. Both must take a standard form:
$r_1$ must accept arguments `w₀`, `p` (parameters not used in this problem), and `t`; $r_2$ must accept arguments `p` and `t`. We will implement these in in-place form to make
it more efficient. $r_1$ will return rate-of-change data of the same type as `w₀`
and $r_2$ will return data `dz` of the same type as `z`

```@repl march
diffusion_rhs!(dw::Nodes,w::Nodes,p,t) = fill!(dw,0.0) # this is r1
boundary_constraint_rhs!(dz::ScalarData,p,t) = fill!(dz,1.0) # this is r2, and sets uniformly to 1
```

Construct the regularization and interpolation operators in their usual
symmetric form, and then set up routines that will provide these operators inside the integrator:

```@repl march
reg = Regularize(X,Δx;issymmetric=true)
Hmat, Emat = RegularizationMatrix(reg,z,w₀)

boundary_constraint_force!(dw::Nodes,z::ScalarData,p) = dw .= Hmat*z # This is B1T
boundary_constraint_op!(dz::ScalarData,y::Nodes,p) = dz .= Emat*y;  # This is B2
```

Note that these last two functions are also in-place, and return data of the same
respective types as $r_1$ and $r_2$.

All of these are assembled into a single `ConstrainedODEFunction`:

```@repl march
f = ConstrainedODEFunction(diffusion_rhs!,boundary_constraint_rhs!,boundary_constraint_force!,
                          boundary_constraint_op!,L,
                          _func_cache=u₀)
```

With the last argument, we supplied a cache variable to enable evaluation of this function.

Now set up the problem, using the same basic notation as in [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl).

```@repl march
tspan = (0.0,20.0)
prob = ODEProblem(f,u₀,tspan)

```

Now solve it. We will set the time-step size to a large value ($1.0$) for demonstration purposes. The method remains stable for any choice.
```@repl march
Δt = 1.0
sol = solve(prob,IFHEEuler(),dt=Δt);
```

Now let's plot it

```@repl march
xg, yg = coordinates(w₀,dx=Δx);
plot(xg,yg,state(sol.u[end]))
plot!(xb,yb,linecolor=:black,linewidth=1.5)
```
![](ifherk.svg)

From a side view, we can see that it enforces the boundary condition:

```@repl march
plot(xg,state(sol.u[end])[65,:],xlabel="x",ylabel="u(x,1)")
savefig("ifherk-side.svg"); nothing # hide
```
![](ifherk-side.svg)
