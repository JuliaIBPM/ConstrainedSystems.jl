
using ConstrainedSystems
using Test
##using TestSetExtensions


#@test isempty(detect_ambiguities(ViscousFlow))
include("saddle.jl")
include("timemarching.jl")


#@testset ExtendedTestSet "All tests" begin
#    @includetests ARGS
#end

#if isempty(ARGS)
#    include("../docs/make.jl")
#end
