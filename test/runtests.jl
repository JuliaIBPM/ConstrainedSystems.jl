
using ConstrainedSystems
using Test
using CartesianGrids
using LinearAlgebra
using Literate

import ConstrainedSystems: recursivecopy!, needs_iteration, ArrayPartition


const GROUP = get(ENV, "GROUP", "All")

ENV["GKSwstype"] = "nul" # removes GKS warnings during plotting

macro mysafetestset(args...)
    name, expr = args
    quote
        ex = quote
          name_str = $$(QuoteNode(name))
          expr_str = $$(QuoteNode(expr))
          mod = gensym(name_str)
          ex2 = quote
              @eval module $mod
                      using Test
                      @testset $name_str $expr_str
                    end
              nothing
          end
          eval(ex2)
        end
        eval(ex)
    end
end

notebookdir = "../examples"
docdir = "../docs/src/manual"
litdir = "./literate"



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

if GROUP == "All" || GROUP == "Aqua"
  include("aqua.jl")
end

if GROUP == "All" || GROUP == "Literate"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      global file_str = "$file"
      global body = :(begin include(joinpath($root,$file)) end)
      #endswith(file,".jl") && startswith(file,"saddle") && @mysafetestset file_str body
      endswith(file,".jl") && @mysafetestset file_str body
    end
  end
end

if GROUP == "Notebooks"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      #endswith(file,".jl") && startswith(file,"saddle") && Literate.notebook(joinpath(root, file),notebookdir)
      endswith(file,".jl") && Literate.notebook(joinpath(root, file),notebookdir)
    end
  end
end

if GROUP == "Documentation"
  for (root, dirs, files) in walkdir(litdir)
    for file in files
      endswith(file,".jl") && Literate.markdown(joinpath(root, file),docdir)
    end
  end
end
