var documenterSearchIndex = {"docs":
[{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"EditURL = \"../../../test/literate/saddlesystems.jl\"","category":"page"},{"location":"manual/saddlesystems/#Saddle-point-systems","page":"Saddle point systems","title":"Saddle point systems","text":"","category":"section"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"CurrentModule = ConstrainedSystems","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"defddt1fracmathrmd1mathrmdt\nrenewcommandvecboldsymbol\nnewcommanduvec1vechat1\nnewcommandutangentuvectau\nnewcommandunormaluvecn\nrenewcommanddmathrmd","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Saddle systems comprise an important part of solving mechanics problems with constraints. In such problems, there is an underlying system to solve, and the addition of constraints requires that the system is subjected to additional forces (constraint forces, or Lagrange multipliers) that enforce these constraints in the system. Examples of such constrained systems are the divergence-free velocity constraint in incompressible flow (for which pressure is the associated Lagrange multiplier field), the no-slip and/or no-flow-through condition in general fluid systems adjacent to impenetrable bodies, and joint constraints in rigid-body mechanics.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"A general saddle-point system has the form","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"left beginarraycc A  B_1^T  B_2  Cendarrayright left(beginarraycuf endarrayright) = left(beginarraycr_1r_2 endarrayright)","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"We are primarily interested in cases when the operator A is symmetric and positive semi-definite, which is fairly typical. It is also fairly common for B_1 = B_2, so that the whole system is symmetric.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"ConstrainedSystems.jl allows us to solve such systems for u and f in a fairly easy way. We need only to provide rules for how to evaluate the actions of the various operators in the system. Let us use an example to show how this can be done.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"using ConstrainedSystems\nusing CartesianGrids\nusing Plots","category":"page"},{"location":"manual/saddlesystems/#Translating-cylinder-in-potential-flow","page":"Saddle point systems","title":"Translating cylinder in potential flow","text":"","category":"section"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"In irrotational, incompressible flow, the streamfunction psi satisfies Laplace's equation,","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"nabla^2 psi = 0","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"On the surface of an impenetrable body, the streamfunction must obey the constraint","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"psi = psi_b","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"where psi_b is the streamfunction associated with the body's motion. Let us suppose the body is moving vertically with velocity 1. Then psi_b = -x for all points inside or on the surface of the body. Thus, the streamfunction field outside this body is governed by Laplace's equation subject to the constraint.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Let us solve this problem on a staggered grid, using the tools discussed in CartesianGrids, including the regularization and interpolation methods to immerse the body shape on the grid. Then our saddle-point system has the form","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"left beginarraycc L  R  E  0endarrayright left(beginarraycpsif endarrayright) = left(beginarrayc0psi_b endarrayright)","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"where L is the discrete Laplacian, R is the regularization operator, and E is the interpolation operator.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Physically, f isn't really a force here, but rather, represents the strengths of distributed singularities on the surface. In fact, this strength represents the jump in normal derivative of psi across the surface. Since this normal derivative is equivalent to the tangential velocity, f is the strength of the bound vortex sheet on the surface. This will be useful to know when we check the value of f obtained in our solution.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"First, let us set up the body, centered at (11) and of radius 12. We will also initialize a data structure for the force:","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"n = 128; θ = range(0,stop=2π,length=n+1);\nxb = 1.0 .+ 0.5*cos.(θ[1:n]); yb = 1.0 .+ 0.5*sin.(θ[1:n]);\nX = VectorData(xb,yb);\nψb = ScalarData(X);\nf = similar(ψb);\nnothing #hide","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Now let's set up a grid of size 102times 102 (including the usual layer of ghost cells) and physical dimensions 2times 2.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"nx = 102; ny = 102; Lx = 2.0; dx = Lx/(nx-2);\nw = Nodes(Dual,(nx,ny));\nψ = similar(w);\nnothing #hide","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"We need to set up the operators now. First, the Laplacian:","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"L = plan_laplacian(size(w),with_inverse=true)","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Note that we have made sure that this operator has an inverse. It is important that this operator, which represents the A matrix in our saddle system, comes with an associated backslash \\\\ operation to carry out the inverse.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Now we need to set up the regularization R and interpolation E operators.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"regop = Regularize(X,dx;issymmetric=true)\nRmat, Emat = RegularizationMatrix(regop,ψb,w);\nnothing #hide","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Now we are ready to set up the system. The solution and right-hand side vectors are set up using SaddleVector.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"rhs = SaddleVector(w,ψb)\nsol = SaddleVector(ψ,f)","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"and the saddle system is then set up with the three operators; the C operator is presumed to be zero when it is not provided.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"A = SaddleSystem(L,Emat,Rmat,rhs)","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Note that all of the operators we have provided are either matrices (like Emat and Rmat) or functions or function-like operators (like L). The SaddleSystem constructor allows either. However, the order is important: we must supply A, B_2, B_1^T, and possibly C, in that order.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Let's solve the system. We need to set the right-hand side. We will set ψb, but this will also change rhs, since that vector is pointing to the same object.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"ψb .= -(xb.-1);\nnothing #hide","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"The right-hand side of the Laplace equation is zero. The right-hand side of the constraint is the specified streamfunction on the body. Note that we have subtracted the circle center from the x positions on the body. The reason for this will be discussed in a moment.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"We solve the system with the convenient shorthand of the backslash:","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"sol .= A\\rhs # hide\n@time sol .= A\\rhs","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Just to point out how fast it can be, we have also timed it. It's pretty fast.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"We can obtain the state vector and the constraint vector from sol using some convenience functions state(sol) and constraint(sol).","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"Now, let's plot the solution in physical space. We'll plot the body shape for reference, also.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"xg, yg = coordinates(w,dx=dx)\nplot(xg,yg,state(sol),xlim=(-Inf,Inf),ylim=(-Inf,Inf))\nplot!(xb,yb,fillcolor=:black,fillrange=0,fillalpha=0.25,linecolor=:black)","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"The solution shows the streamlines for a circle in vertical motion, as expected. All of the streamlines inside the circle are vertical.","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"","category":"page"},{"location":"manual/saddlesystems/","page":"Saddle point systems","title":"Saddle point systems","text":"This page was generated using Literate.jl.","category":"page"},{"location":"manual/methods/#Index","page":"Index","title":"Index","text":"","category":"section"},{"location":"manual/methods/","page":"Index","title":"Index","text":"DocTestSetup = quote\n  using ConstrainedSystems\nend","category":"page"},{"location":"manual/methods/","page":"Index","title":"Index","text":"Modules  = [ConstrainedSystems]\nOrder   = [:type, :function]","category":"page"},{"location":"manual/methods/#ConstrainedSystems.ConstrainedODEFunction","page":"Index","title":"ConstrainedSystems.ConstrainedODEFunction","text":"    ConstrainedODEFunction(r1,r2,B1,B2[,L][,C])\n\nThis specifies the functions and operators that comprise an ODE problem with the form\n\ndfracdydt = Ly - B_1 z + r_1(yt)\n\nB_2 y + C z = r_2(xt)\n\nwhere y is the state, z is a constraint force, and x is an auxiliary state describing the constraints.\n\nThe optional linear operator L defaults to zeros. The B1 and B2 functions must be of the respective in-place forms B1(dy,z,x,p) (to compute the action of B1 on z) and B2(dz,y,x,p) (to compute the action of B2 on y). The function r1 must of the in-place form r1(dy,y,x,p,t), and r2 must be in the in-place form r2(dz,x,p,t). The C function can be omitted, but if it is included, then it must be of the form C(dz,z,x,p) (to compute the action of C on z). Alternatively, one can supply out-of-place forms, respectively, as B1(z,x,p), B2(y,x,p), C(z,x,p), r1(y,x,p,t) and r2(x,p,t).\n\nAn optional keyword argument param_update_func can be used to set a function that updates problem parameters with the current solution. This function must take the in-place form f(q,u,p,t) or out of place form f(u,p,t) to create some q based on u, where y = state(u), z = constraint(u) and x = aux_state(u). (Note that q might enter the function simply as p, to be mutated.) This function can be used to update B1, B2, and C, for example.\n\nWe can also include another (unconstrained) set of equations to the set above in order to update x:\n\ndfracdxdt = r_1x(upt)\n\nIn this case, the right-hand side has access to the entire u vector. We would pass the pair of r1 functions as an ArrayPartition.\n\n\n\n\n\n","category":"type"},{"location":"manual/methods/#ConstrainedSystems.SaddleSystem-Union{Tuple{TS}, Tuple{T}, Tuple{LinearMaps.LinearMap{T}, LinearMaps.LinearMap{T}, LinearMaps.LinearMap{T}, LinearMaps.LinearMap{T}, LinearMaps.LinearMap{T}, LinearMaps.LinearMap{T}, Any, Any}} where {T, TS<:SchurSolverType}","page":"Index","title":"ConstrainedSystems.SaddleSystem","text":"SaddleSystem\n\nConstruct a saddle-point system operator from the constituent operator blocks. The resulting object can be used with * and \\ to multiply and solve. The saddle-point problem has the form\n\nbeginbmatrixA  B_1^T  B_2  C endbmatrix beginpmatrix u  f endpmatrix = beginpmatrix r_1  r_2 endpmatrix\n\nConstructors\n\nSaddleSystem(A::AbstractMatrix,B₂::AbstractMatrix,B₁ᵀ::AbstractMatrix,C::AbstractMatrix[,eltype=Float64]). Blocks are given as matrices. Must have consistent sizes to stack appropriately. If this is called with SaddleSystem(A,B₂,B₁ᵀ), it sets C to zero automatically.\n\nSaddleSystem(A,B₂,B₁ᵀ,C,u,f[,eltype=Float64]). Operators A, B₂, B₁ᵀ, C are given in various forms, including matrices, functions, and function-like objects. u and f are examples of the data types in the corresponding solution and right-hand side vectors. Guidelines:\n\nThe entries A and B₂ must be able to act upon u (either by multiplication or as a function) and B₁ᵀ and C must be able to act on f (also, either by multiplication or as a function).\nA and B₁ᵀ should return data of type u, and B₂ and C should return data of type f.\nA must be invertible and be outfitted with operators `andldiv!`.\nBoth u and f must be subtypes of AbstractArray: they must be equipped with size and vec functions and with a constructor of the form T(data) where T is the data type of u or f and data is the wrapped data array.\n\nIf called as SaddleSystem(A,B₂,B₁ᵀ,u,f), the C block is omitted and assumed to be zero.\n\nIf called with SaddleSystem(A,u), this is equivalent to calling SaddleSystem(A,nothing,nothing,u,[]), then this reverts to the unconstrained system described by operator A.\n\nThe list of vectors u and f in any of these constructors can be bundled together as a SaddleVector, e.g. SaddleSystem(A,B₂,B₁ᵀ,SaddleVector(u,f)).\n\nAn optional keyword argument solver= can be used to specify the type of solution for the Schur complement system. By default, this is set to Direct, and the Schur complement matrix is formed, factorized, and stored. This can be changed to Iterative, in which case an iterative solver is determined by the linsolve function of KrylovKit.jl.\n\n\n\n\n\n","category":"method"},{"location":"manual/methods/#ConstrainedSystems.SaddleVector","page":"Index","title":"ConstrainedSystems.SaddleVector","text":"SaddleVector(u,f)\n\nConstruct a vector of a state part u and constraint part f of a saddle-point vector, to be associated with a SaddleSystem.\n\n\n\n\n\n","category":"type"},{"location":"manual/methods/#Base.eltype-Union{Tuple{SaddleSystem{T, Ns, Nc, TU, TF, TS} where {TU, TF, TS<:SchurSolverType}}, Tuple{Nc}, Tuple{Ns}, Tuple{T}} where {T, Ns, Nc}","page":"Index","title":"Base.eltype","text":"Base.eltype(::SaddleSystem)\n\nReport the element type of a SaddleSystem.\n\n\n\n\n\n","category":"method"},{"location":"manual/methods/#Base.size-Union{Tuple{SaddleSystem{T, Ns, Nc, TU, TF, TS} where {TU, TF, TS<:SchurSolverType}}, Tuple{Nc}, Tuple{Ns}, Tuple{T}} where {T, Ns, Nc}","page":"Index","title":"Base.size","text":"Base.size(::SaddleSystem)\n\nReport the size of a SaddleSystem.\n\n\n\n\n\n","category":"method"},{"location":"manual/methods/#ConstrainedSystems.aux_state-Tuple{Any}","page":"Index","title":"ConstrainedSystems.aux_state","text":"aux_state(x)\n\nProvide the auxiliary state part of the given vector x\n\n\n\n\n\n","category":"method"},{"location":"manual/methods/#ConstrainedSystems.constraint-Tuple{RecursiveArrayTools.ArrayPartition}","page":"Index","title":"ConstrainedSystems.constraint","text":"constraint(x::SaddleVector)\n\nProvide the constraint part of the given saddle vector x\n\n\n\n\n\n","category":"method"},{"location":"manual/methods/#ConstrainedSystems.r1vector-Tuple{}","page":"Index","title":"ConstrainedSystems.r1vector","text":"r1vector([;state_r1=][,aux_r1=])\n\nBuild a vector of the r1 functions for the state ODEs and auxiliary state ODEs.\n\n\n\n\n\n","category":"method"},{"location":"manual/methods/#ConstrainedSystems.solvector-Tuple{}","page":"Index","title":"ConstrainedSystems.solvector","text":"solvector([;state=][,constraint=][,aux_state=])\n\nBuild a solution vector for a constrained system. This takes three optional keyword arguments: state, constraint, and aux_state. If only a state is supplied, then the constraint is set to an empty vector and the system is assumed to correspond to an unconstrained system. (aux_state is ignored in this situation.)\n\n\n\n\n\n","category":"method"},{"location":"manual/methods/#ConstrainedSystems.state-Tuple{RecursiveArrayTools.ArrayPartition}","page":"Index","title":"ConstrainedSystems.state","text":"state(x::SaddleVector)\n\nProvide the state part of the given saddle vector x\n\n\n\n\n\n","category":"method"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"EditURL = \"../../../test/literate/timemarching.jl\"","category":"page"},{"location":"manual/timemarching/#Time-marching","page":"Time marching","title":"Time marching","text":"","category":"section"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"CurrentModule = ConstrainedSystems","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"defddt1fracmathrmd1mathrmdt\nrenewcommandvecboldsymbol\nnewcommanduvec1vechat1\nnewcommandutangentuvectau\nnewcommandunormaluvecn\nrenewcommanddmathrmd","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"ConstrainedSystems.jl is equipped with tools for solving systems of equations of the general form of half-explicit differential-algebraic equations,","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"ddt y = L y - B_1^T(yt) z + r_1(yt) quad B_2(yt) y + C(yt) z = r_2(t) quad y(0) = y_0","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"where z is the Lagrange multiplier for enforcing the constraints on y. Note that the constraint operators may depend on the state and on time. The linear operator L may be a matrix or a scalar, but is generally independent of time. (The method of integrating factors can deal with time-dependent L, but we don't encounter such systems in the constrained systems context so we won't discuss them.) Our objective is to solve for y(t) and z(t).","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"using ConstrainedSystems\nusing CartesianGrids\nusing Plots","category":"page"},{"location":"manual/timemarching/#Constrained-integrating-factor-systems","page":"Time marching","title":"Constrained integrating factor systems","text":"","category":"section"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Let's demonstrate this on the example of heat diffusion from a circular ring whose temperature is held constant. In this case, L is the discrete Laplace operator times the heat diffusivity, r_1 is zero (in the absence of volumetric heating sources), and r_2 is the temperature of the ring. The operators B_1^T and B_2 will be the regularization and interpolation operators between discrete point-wise data on the ring and the field data. We will also include a C operator that slightly regularizes the constraint.","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"The ring will have radius 12 and fixed temperature 1, and the heat diffusivity is 1. (In other words, the problem has been non-dimensionalized by the diameter of the circle, the dimensional ring temperature, and the dimensional diffusivity.)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"First, we will construct a field to accept the temperature on","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"nx = 129; ny = 129; Lx = 2.0; Δx = Lx/(nx-2);\nw₀ = Nodes(Dual,(nx,ny)); # field initial condition\nnothing #hide","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now set up a ring of points on the circle at center (11).","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"n = 128; θ = range(0,stop=2π,length=n+1);\nR = 0.5; xb = 1.0 .+ R*cos.(θ); yb = 1.0 .+ R*sin.(θ);\nX = VectorData(xb[1:n],yb[1:n]);\nz = ScalarData(X); # to be used as the Lagrange multiplier\nnothing #hide","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Together, w₀ and z comprise the initial solution vector:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"u₀ = solvector(state=w₀,constraint=z)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now set up the operators. We first set up the linear operator, a Laplacian endowed with its inverse:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"L = plan_laplacian(w₀,with_inverse=true)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now the right-hand side operators for the ODEs and constraints. Both must take a standard form: r_1 must accept arguments w₀, p (parameters not used in this problem), and t; r_2 must accept arguments p and t. We will implement these in in-place form to make it more efficient. r_1 will return rate-of-change data of the same type as w₀ and r_2 will return data dz of the same type as z","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"diffusion_rhs!(dw::Nodes,w::Nodes,x,p,t) = fill!(dw,0.0) # this is r1\nboundary_constraint_rhs!(dz::ScalarData,x,p,t) = fill!(dz,1.0) # this is r2, and sets uniformly to 1","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Construct the regularization and interpolation operators in their usual symmetric form, and then set up routines that will provide these operators inside the integrator:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"reg = Regularize(X,Δx;issymmetric=true)\nHmat, Emat = RegularizationMatrix(reg,z,w₀)\n\nboundary_constraint_force!(dw::Nodes,z::ScalarData,x,p) = dw .= Hmat*z # This is B1T\nboundary_constraint_op!(dz::ScalarData,y::Nodes,x,p) = dz .= Emat*y;  # This is B2\nnothing #hide","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Construct a constraint regularization operator (the C operator)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"boundary_constraint_reg!(dz::ScalarData,z::ScalarData,x,p) = dz .= -0.1*z; # This is C\nnothing #hide","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Note that these last two functions are also in-place, and return data of the same respective types as r_1 and r_2.","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"All of these are assembled into a single ConstrainedODEFunction:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"f = ConstrainedODEFunction(diffusion_rhs!,boundary_constraint_rhs!,boundary_constraint_force!,\n                          boundary_constraint_op!,L,boundary_constraint_reg!,_func_cache=u₀)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"With the last argument, we supplied a cache variable to enable evaluation of this function.","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now set up the problem, using the same basic notation as in DifferentialEquations.jl.","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"tspan = (0.0,20.0)\nprob = ODEProblem(f,u₀,tspan)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now solve it. We will set the time-step size to a large value (10) for demonstration purposes. The method remains stable for any choice.","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Δt = 1.0\nsol = solve(prob,IFHEEuler(),dt=Δt);\nnothing #hide","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now let's plot it","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"xg, yg = coordinates(w₀,dx=Δx);\nplot(xg,yg,state(sol.u[end]))\nplot!(xb,yb,linecolor=:black,linewidth=1.5)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"From a side view, we can see that it enforces the boundary condition:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"plot(xg,state(sol.u[end])[65,:],xlabel=\"x\",ylabel=\"u(x,1)\")","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"The Lagrange multiplier distribution is nearly uniform","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"plot(constraint(sol.u[end]),ylim=(-0.5,0))","category":"page"},{"location":"manual/timemarching/#Systems-with-variable-constraints","page":"Time marching","title":"Systems with variable constraints","text":"","category":"section"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"In some cases, the constraint operators may vary with the state vector. A good example of this is a swinging pendulum, with its equations expressed in Cartesian coordinates. The constraint we wish to enforce is that the length of the pendulum is constant: x^2+y^2 = l^2. Though not mathematically necessary, it also helps to enforce a tangency condition, xu + yv = 0, where u and v are the rates of change of x and y. Note that this is simply the derivative of the first constraint. (If expressed in cylindrical coordinates, the constraint is enforced automatically, simply by expressing the equations for theta.)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"The governing equations are","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"ddt x = u - x mu  ddt y = v - y mu   ddt u = -x lambda   ddt v = -g - ylambda","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"with Lagrange multipliers mu and lambda, and the constraints are","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"x^2+y^2 = l^2 xu + yu = 0","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"These are equivalently expressed as","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"left beginarraycccc x  y  0  0endarrayrightleft beginarrayc x  y  u  v endarrayright = l^2","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"and","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"left beginarraycccc 0  0  x  yendarrayrightleft beginarrayc x  y  u  v endarrayright = 0","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"The operators B_1^T and B_2 are thus","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"B_1^T = left beginarraycc x  0  y  0  0  x  0  y endarrayright","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"and","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"B_2 = left beginarraycccc x  y  0  0  0  0  x  y endarrayright","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"That is, the operators are dependent on the state. In this package, we handle this by providing a parameter that can be dynamically updated. We will get to that later. First, let's set up the physical parameters","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"l = 1.0\ng = 1.0\nparams = [l,g]","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"and initial condition:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"θ₀ = π/2\ny₀ = Float64[l*sin(θ₀),-l*cos(θ₀),0,0]\nz₀ = Float64[0.0, 0.0] # Lagrange multipliers\nu₀ = solvector(state=y₀,constraint=z₀)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now, we will set up the basic form of the constraint operators and assemble these with the other parameters with the help of a type we'll define here:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"struct ProblemParams{P,BT1,BT2}\n    params :: P\n    B₁ᵀ :: BT1\n    B₂ :: BT2\nend\n\nB1T = zeros(4,2) # set to zeros for now\nB2 = zeros(2,4)  # set to zeros for now\np₀ = ProblemParams(params,B1T,B2);\nnothing #hide","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"We will now define the operators of the problem, all in in-place form:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"function pendulum_rhs!(dy::Vector{Float64},y::Vector{Float64},x,p,t)\n    dy[1] = y[3]\n    dy[2] = y[4]\n    dy[3] = 0.0\n    dy[4] = -p.params[2]\n    return dy\nend # r1\n\nfunction length_constraint_rhs!(dz::Vector{Float64},x,p,t)\n    dz[1] = p.params[1]^2\n    dz[2] = 0.0\n    return dz\nend # r2","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"The B1 function. This returns B1*z. It uses an existing B1 supplied by p.","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"function length_constraint_force!(dy::Vector{Float64},z::Vector{Float64},x,p)\n    dy .= p.B₁ᵀ*z\nend","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"The B2 function. This returns B2*y. It uses an existing B2 supplied by p.","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"function length_constraint_op!(dz::Vector{Float64},y::Vector{Float64},x,p)\n    dz .= p.B₂*y\nend","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now, we need to provide a means of updating the parameter structure with the current state of the system. This is done in-place, just as for the other operators:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"function update_p!(q,u,p,t)\n    y = state(u)\n    fill!(q.B₁ᵀ,0.0)\n    fill!(q.B₂,0.0)\n    q.B₁ᵀ[1,1] = y[1]; q.B₁ᵀ[2,1] = y[2]; q.B₁ᵀ[3,2] = y[1]; q.B₁ᵀ[4,2] = y[2]\n    q.B₂[1,1] = y[1]; q.B₂[1,2] = y[2]; q.B₂[2,3] = y[1]; q.B₂[2,4] = y[2]\n    return q\nend","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Finally, assemble all of them together:","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"f = ConstrainedODEFunction(pendulum_rhs!,length_constraint_rhs!,length_constraint_force!,\n                                length_constraint_op!,\n                               _func_cache=deepcopy(u₀),param_update_func=update_p!)","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Now solve the system","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"tspan = (0.0,10.0)\nprob = ODEProblem(f,u₀,tspan,p₀)\n\nΔt = 1e-2\nsol = solve(prob,LiskaIFHERK(),dt=Δt);\nnothing #hide","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"Plot the solution","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"plot(sol.t,sol[1,:],label=\"x\",xlabel=\"t\")\nplot!(sol.t,sol[2,:],label=\"y\")","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"and here is the trajectory","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"plot(sol[1,:],sol[2,:],ratio=1,legend=:false,title=\"Trajectory\",xlabel=\"x\",ylabel=\"y\")","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"","category":"page"},{"location":"manual/timemarching/","page":"Time marching","title":"Time marching","text":"This page was generated using Literate.jl.","category":"page"},{"location":"#ConstrainedSystems.jl","page":"Home","title":"ConstrainedSystems.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"tools for solving constrained dynamical systems","category":"page"},{"location":"","page":"Home","title":"Home","text":"defddt1fracmathrmd1mathrmdt\n\nrenewcommandvecboldsymbol\nnewcommanduvec1vechat1\nnewcommandutangentuvectau\nnewcommandunormaluvecn\n\nrenewcommanddmathrmd","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package contains several tools for solving and advancing (large-scale) dynamical systems with constraints. These systems generically have the form","category":"page"},{"location":"","page":"Home","title":"Home","text":"ddty = L u - B_1^T z + r_1(yt) quad B_2 y = r_2(t) quad y(0) = y_0","category":"page"},{"location":"","page":"Home","title":"Home","text":"where y is a state vector, L is a linear operator with an associated matrix exponential (integrating factor), and z is a constraint force vector (i.e., Lagrange multipliers). Systems of this type might arise from, e.g., incompressible fluid dynamics, rigid-body mechanics, or couplings of such systems.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Some of the key components of this package are","category":"page"},{"location":"","page":"Home","title":"Home","text":"Tools for solving linear algebra problems with constraints and associated Lagrange multipliers, known generically as saddle point systems. The sizes of these systems might be large.\nTime integrators that can incorporate these constraints, such as half-explicit Runge-Kutta (HERK) and integrating factor Runge-Kutta (IFRK), or their combination (IF-HERK). These","category":"page"},{"location":"","page":"Home","title":"Home","text":"extend the tools in the DifferentialEquations.jl package, and utilize the same basic syntax for setting up a problem and solving it.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Allowance for variable constraint operators B_1^T and B_2,","category":"page"},{"location":"","page":"Home","title":"Home","text":"through the use of a variable parameter argument and an associated parameter update function.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The ability to add an auxiliary (unconstrained) system of equations and state that the","category":"page"},{"location":"","page":"Home","title":"Home","text":"constraint operators B_1^T and B_2 and the right-hand side r_2 depend upon.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package works on Julia 1.6 and above and is registered in the general Julia registry. To install from the REPL, type e.g.,","category":"page"},{"location":"","page":"Home","title":"Home","text":"] add ConstrainedSystems","category":"page"},{"location":"","page":"Home","title":"Home","text":"Then, in any version, type","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using ConstrainedSystems","category":"page"},{"location":"","page":"Home","title":"Home","text":"The plots in this documentation are generated using Plots.jl. You might want to install that, too, to follow the examples.","category":"page"}]
}
