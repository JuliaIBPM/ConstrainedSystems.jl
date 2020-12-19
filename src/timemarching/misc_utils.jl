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

function recursivecopy!(dest :: T, src :: T) where {T}
    fields = fieldnames(T)
    for f in fields
        tmp = getfield(dest,f)
        copy!(tmp,getfield(src,f))
    end
    dest
end

function recursivecopy!(dest :: AbstractArray{T}, src :: AbstractArray{T}) where {T}
    for i in eachindex(src)
        recursivecopy!(dest[i],src[i])
    end
    return dest
end
