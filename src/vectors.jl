### Right-hand side and solution vectors

const SaddleVector = ArrayPartition

"""
    SaddleVector(u,f)

Construct a vector of a state part `u` and constraint part `f` of a
saddle-point vector, to be associated with a [`SaddleSystem`](@ref).
"""
function SaddleVector end


"""
    solvector([;state=][,constraint=][,aux_state=])

Build a solution vector for a constrained system. This takes three optional keyword
arguments: `state`, `constraint`, and `aux_state`. If only a state is supplied,
then the constraint is set to an empty vector and the system is assumed to
correspond to an unconstrained system. (`aux_state` is ignored in this situation.)
"""
solvector(;state=nothing,constraint=nothing,aux_state=nothing) = _solvector(state,constraint,aux_state)
_solvector(::Nothing,::Nothing,::Nothing) = nothing
_solvector(s,::Nothing,::Nothing) = ArrayPartition(s,empty(s))
_solvector(s,::Nothing,aux) = ArrayPartition(s,empty(s))
_solvector(s,c,::Nothing) = ArrayPartition(s,c)
#_solvector(s,c,aux) = ArrayPartition(_solvector(s,c,nothing),aux)
_solvector(s,c,aux) = ArrayPartition(s,c,aux)


mainvector(u) = u
mainvector(u::ArrayPartition) = ArrayPartition(u.x[1],u.x[2])




"""
    state(x::SaddleVector)

Provide the state part of the given saddle vector `x`
"""
state(u::ArrayPartition) = mainvector(u).x[1]
state(u) = u

"""
    constraint(x::SaddleVector)

Provide the constraint part of the given saddle vector `x`
"""
constraint(u::ArrayPartition) = mainvector(u).x[2]
constraint(u) = eltype(u)[]


"""
    aux_state(x)

Provide the auxiliary state part of the given vector `x`
"""
aux_state(u) = nothing # Array{eltype(u)}(undef,0,0)
aux_state(u::ArrayPartition{T,Tuple{F1,F2,F3}}) where {T,F1,F2,F3} = u.x[3]



for f in (:state,:constraint,:aux_state)
  @eval $f(a::AbstractArray{T}) where {T<:ArrayPartition} = map($f,a)
end


#SaddleVector(u::TU,f::TF) where {TU,TF} = ArrayPartition(u,f)




"""
    r1vector([;state_r1=][,aux_r1=])

Build a vector of the `r1` functions for the state ODEs and auxiliary state ODEs.
"""
r1vector(;state_r1=nothing,aux_r1=nothing) = _r1vector(state_r1,aux_r1)
_r1vector(::Nothing,::Nothing) = nothing
_r1vector(s,::Nothing) = s
_r1vector(::Nothing,a) = nothing
_r1vector(s,a) = ArrayPartition((s,a))

hasaux(r1) = false
hasaux(r1::ArrayPartition) = true

state_r1(r1) = r1
state_r1(r1::ArrayPartition) = r1.x[1]

aux_r1(r1) = nothing
aux_r1(r1::ArrayPartition) = r1.x[2]
