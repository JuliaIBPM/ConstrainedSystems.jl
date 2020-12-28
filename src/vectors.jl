### Right-hand side and solution vectors

const SaddleVector = ArrayPartition

"""
    solvector()

Build a solution vector for a constrained system. This takes three optional keyword
arguments: `state`, `constraint`, and `aux_state`.
"""
solvector(;state=nothing,constraint=nothing,aux_state=nothing) = _solvector(state,constraint,aux_state)
_solvector(::Nothing,::Nothing,::Nothing) = nothing
_solvector(s,::Nothing,::Nothing) = s
_solvector(s,::Nothing,aux) = s
_solvector(s,c,::Nothing) = ArrayPartition(s,c)
#_solvector(s,c,aux) = ArrayPartition(_solvector(s,c,nothing),aux)
_solvector(s,c,aux) = ArrayPartition(s,c,aux)


mainvector(u) = u
#mainvector(u::ArrayPartition{T,Tuple{A,F}}) where {T,A<:ArrayPartition,F} = u.x[1]
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
#aux_state(u::ArrayPartition{T,Tuple{A,F}}) where {T,A<:ArrayPartition,F} = u.x[2]
aux_state(u::ArrayPartition{T,NTuple{3,F}}) where {T,F} = u.x[3]



for f in (:state,:constraint,:aux_state)
  @eval $f(a::AbstractArray{T}) where {T<:ArrayPartition} = map($f,a)
end


#SaddleVector(u::TU,f::TF) where {TU,TF} = ArrayPartition(u,f)

"""
    SaddleVector(u,f)

Construct a vector of a state part `u` and constraint part `f` of a
saddle-point vector, to be associated with a [`SaddleSystem`](@ref).
"""
function SaddleVector end
