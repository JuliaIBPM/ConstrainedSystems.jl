#=
This is pirated from
https://github.com/SciML/OrdinaryDiffEq.jl/blob/0e38705ac4ace80a51961d689ff64bea6b3bce73/src/misc_utils.jl#L20
because for some reason it is not working by importing it
=#
macro cache(expr)
  name = expr.args[2].args[1].args[1]
  fields = [x for x in expr.args[3].args if typeof(x)!=LineNumberNode]
  cache_vars = Expr[]
  jac_vars = Pair{Symbol,Expr}[]
  for x in fields
    if x.args[2] == :uType || x.args[2] == :rateType ||
       x.args[2] == :kType || x.args[2] == :uNoUnitsType
      push!(cache_vars,:(c.$(x.args[1])))
    elseif x.args[2] == :DiffCacheType
      push!(cache_vars,:(c.$(x.args[1]).du))
      push!(cache_vars,:(c.$(x.args[1]).dual_du))
    end
  end
  quote
    $expr
    $(esc(:full_cache))(c::$name) = tuple($(cache_vars...))
  end
end

# Need to allow UNITLESS_ABS2 to work on empty vectors
# Should add this into DiffEqBase
#@inline UNITLESS_ABS2(x::AbstractArray) = (isempty(x) && return sum(UNITLESS_ABS2,zero(eltype(x))); sum(UNITLESS_ABS2, x))

# A workaround that avoids redefinition
import OrdinaryDiffEq.DiffEqBase: UNITLESS_ABS2, ODE_DEFAULT_NORM, recursive_length
@inline MY_UNITLESS_ABS2(x::Number) = abs2(x)
@inline MY_UNITLESS_ABS2(x::AbstractArray) = (isempty(x) && return sum(MY_UNITLESS_ABS2,zero(eltype(x))); sum(MY_UNITLESS_ABS2, x))
@inline MY_UNITLESS_ABS2(x::ArrayPartition) = sum(MY_UNITLESS_ABS2, x.x)

@inline ODE_DEFAULT_NORM(u::ArrayPartition,t) = sqrt(MY_UNITLESS_ABS2(u)/recursive_length(u))


@inline _l2norm(u) = sqrt(recursive_mean(map(x -> float(x).^2,u)))


zero_vec!(::Nothing) = nothing
zero_vec!(u) = fill!(u,0.0)
function zero_vec!(u::ArrayPartition)
    for x in u.x
        fill!(x,0.0)
    end
    return nothing
end



function recursivecopy!(dest :: T, src :: T) where {T}
    fields = fieldnames(T)
    for f in fields
        tmp = getfield(dest,f)
        #tmp .= getfield(src,f)
        recursivecopy!(tmp,getfield(src,f))
    end
    dest
end

function recursivecopy!(dest :: AbstractArray{T}, src :: AbstractArray{T}) where {T}
    for i in eachindex(src)
        recursivecopy!(dest[i],src[i])
    end
    return dest
end


# Seed the state vector with two sets of random values, apply the constraint operator on a
needs_iteration(f::ConstrainedODEFunction{iip},u,p,rate_prototype) where {iip} = _needs_iteration(f,u,p,rate_prototype,Val(iip))

function _needs_iteration(f,u,p,rate_prototype,::Val{true})

    pseed = deepcopy(p)
    u_target, useed = (zero(u) for i in 1:2)
    yseed = state(useed)
    y_target = state(u_target)
    fill!(y_target,1.0)

    dutmp, dudiff = (zero(rate_prototype) for i in 1:2)
    dzdiff = constraint(dudiff)

    yseed .= randn(size(yseed))
    f.param_update_func(pseed,useed,p,0.0)
    _constraint_neg_B2!(dutmp,f,u_target,pseed,0.0)
    dudiff .= dutmp

    yseed .= randn(size(yseed))
    f.param_update_func(pseed,useed,p,0.0)
    _constraint_neg_B2!(dutmp,f,u_target,pseed,0.0)
    dudiff .-= dutmp
    !(_l2norm(dzdiff) == 0.0)
end

function _needs_iteration(f,u,p,rate_prototype,::Val{false})

    pseed = deepcopy(p)
    u_target, useed = (zero(u) for i in 1:2)
    yseed = state(useed)
    y_target = state(u_target)
    fill!(y_target,1.0)

    dutmp, dudiff = (zero(rate_prototype) for i in 1:2)
    dzdiff = constraint(dudiff)

    yseed .= randn(size(yseed))
    pseed = f.param_update_func(useed,p,0.0)
    dutmp = _constraint_neg_B2(f,u_target,pseed,0.0)
    dudiff .= dutmp

    yseed .= randn(size(yseed))
    pseed = f.param_update_func(useed,p,0.0)
    dutmp = _constraint_neg_B2(f,u_target,pseed,0.0)
    dudiff .-= dutmp
    !(_l2norm(dzdiff) == 0.0)
end
