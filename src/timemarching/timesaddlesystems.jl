
# called directly by alg_cache for iip algorithms and indirectly by non-sc algorithms
function SaddleSystem(A,f::ConstrainedODEFunction{true},p,pold,ducache,solver;cfact=1.0)
    nully, nullz = state(ducache), constraint(ducache)
    du_aux = aux_state(ducache)
    @inline B₁ᵀ(z) = (zero_vec!(ducache);
                     _ode_neg_B1!(ducache,f,solvector(state=nully,constraint=z,aux_state=du_aux),pold,0.0);
                     state(ducache) .*= -1.0; return state(ducache))
    @inline B₂(y) = (zero_vec!(ducache);
                     _constraint_neg_B2!(ducache,f,solvector(state=y,constraint=nullz,aux_state=du_aux),p,0.0);
                     constraint(ducache) .*= -1.0; return constraint(ducache))
    @inline C(z) = (zero_vec!(ducache);
                    _constraint_neg_C!(ducache,f,solvector(state=nully,constraint=z,aux_state=du_aux),p,0.0);
                     constraint(ducache) .*= -cfact; return constraint(ducache))
    SaddleSystem(A,B₂,B₁ᵀ,C,mainvector(ducache),solver=solver)
end

# called directly by oop algorithms
function SaddleSystem(A,f::ConstrainedODEFunction{false},p,pold,ducache,solver;cfact=1.0)

    nully, nullz = state(ducache), constraint(ducache)
    du_aux = aux_state(ducache)
    @inline B₁ᵀ(z) = (zero_vec!(ducache); ducache .= _ode_neg_B1(f,solvector(state=nully,constraint=z,aux_state=du_aux),pold,0.0);
                     state(ducache) .*= -1.0; return state(ducache))
    @inline B₂(y) = (zero_vec!(ducache);
                     ducache .= _constraint_neg_B2(f,solvector(state=y,constraint=nullz,aux_state=du_aux),p,0.0);
                     constraint(ducache) .*= -1.0; return constraint(ducache))
    @inline C(z) = (zero_vec!(ducache);
                     ducache .= _constraint_neg_C(f,solvector(state=nully,constraint=z,aux_state=du_aux),p,0.0);
                     constraint(ducache) .*= -cfact; return constraint(ducache))
    SaddleSystem(A,B₂,B₁ᵀ,C,mainvector(ducache),solver=solver)
end

# this version is called by in-place algorithms
@inline SaddleSystem(S::SaddleSystem,A,f::ConstrainedODEFunction,p,pold,
                      cache::ConstrainedODEMutableCache{sc,solverType}) where {sc,solverType} =
          SaddleSystem(S,A,f,p,pold,cache.dutmp,solverType,Val(sc))

# non-static constraints
@inline SaddleSystem(S::SaddleSystem,A,f::ConstrainedODEFunction,p,pold,ducache,solver,
                      ::Val{false}) = SaddleSystem(A,f,p,pold,ducache,solver)

# static constraints
@inline SaddleSystem(S::SaddleSystem,A,f::ConstrainedODEFunction,p,pold,ducache,solver,
                      ::Val{true}) = S


function B1_times_z!(u,S::SaddleSystem)
    state(u) .= typeof(state(u))(S.A⁻¹B₁ᵀf)
    return u
end

function B1_times_z(u,S::SaddleSystem)
    out = zero(u)
    state(out) .= typeof(state(u))(S.A⁻¹B₁ᵀf)
    return out
end
