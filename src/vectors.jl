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
_solvector(s,::Nothing,::Nothing) = ArrayPartition(s,_empty(s))
_solvector(s,::Nothing,aux) = ArrayPartition(s,_empty(s))
_solvector(s,c,::Nothing) = ArrayPartition(s,c)
#_solvector(s,c,aux) = ArrayPartition(_solvector(s,c,nothing),aux)
_solvector(s,c,aux) = ArrayPartition(s,c,aux)


mainvector(u) = u
mainvector(u::ArrayPartition) = ArrayPartition(u.x[1],u.x[2])

_empty(s) = Vector{eltype(s)}()


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


# To extend interpolation to ArrayPartition but preserve the types of entries
Base.BroadcastStyle(::Type{<:ArrayPartition}) = Broadcast.ArrayStyle{ArrayPartition}()

# This is necessary for scalar multiplications
function Base.similar(bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{ArrayPartition}}, ::Type{ElType}) where ElType
    ap = find_ap(bc)
    ArrayPartition(similar.(ap.x))
end

# Find the first ArrayPartition in the Broadcast wrapper
find_ap(bc::Base.Broadcast.Broadcasted) = find_ap(bc.args)
find_ap(args::Tuple) = find_ap(find_ap(args[1]), Base.tail(args))
find_ap(x) = x
find_ap(::Tuple{}) = nothing
find_ap(a::ArrayPartition, rest) = a
find_ap(::Any, rest) = find_ap(rest)


# Extended copyto! for broadcasting to ArrayPartition, fusing the dot operations
# It calls copyto! of the individual x fields, since these are possibly optimized
@inline function Base.copyto!(dest::ArrayPartition, bc::Broadcast.Broadcasted{Broadcast.ArrayStyle{ArrayPartition}})
    if bc.f === identity && bc.args isa Tuple{AbstractArray} # only a single input argument to broadcast!
        A = bc.args[1]
        if axes(dest) == axes(A)
            return copyto!(dest, A)
        end
    end
    bcf = Base.Broadcast.flatten(bc)
    for i in eachindex(dest.x)
        copyto!(dest.x[i],unpack(bcf,i))
    end

    return dest
end

# This fixes a flaw in the copyto! in RecursiveArrayTools
function Base.copyto!(dest::ArrayPartition,src::ArrayPartition)
    @assert length(src) == length(dest)
    if size.(dest.x) == size.(src.x)
       @inbounds for i in eachindex(src.x)
        dest.x[i] .= src.x[i]
      end
    else
      cnt = 0
      for i in eachindex(dest.x)
        x = dest.x[i]
        for k in eachindex(x)
          cnt += 1
          x[k] = src[cnt]
        end
      end
    end
    dest
end

# Generate a Broadcasted entry with only one x field entry of ArrayPartition
@inline unpack(bc::Broadcast.Broadcasted, i) = Broadcast.Broadcasted(bc.f, unpack_args(i, bc.args))
unpack(x,::Any) = x
unpack(x::ArrayPartition, i) = x.x[i]

@inline unpack_args(i, args::Tuple) = (unpack(args[1], i), unpack_args(i, Base.tail(args))...)
unpack_args(i, args::Tuple{Any}) = (unpack(args[1], i),)
unpack_args(::Any, args::Tuple{}) = ()
