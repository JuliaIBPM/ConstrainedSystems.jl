
using ConstrainedSystems
using Test
using CartesianGrids
using LinearAlgebra

import ConstrainedSystems: recursivecopy!, needs_iteration, ArrayPartition


const GROUP = get(ENV, "GROUP", "All")

if GROUP == "All" || GROUP == "Utils"
  include("utils.jl")
end

if GROUP == "All" || GROUP == "Types"
  include("types.jl")
end

if GROUP == "All" || GROUP == "Saddle"
  include("saddle.jl")
end

if GROUP == "All" || GROUP == "Convergence"
  include("algconvergence.jl")
end
