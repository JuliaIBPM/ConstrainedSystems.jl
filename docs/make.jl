using Documenter, ConstrainedSystems

makedocs(
    sitename = "ConstrainedSystems.jl",
    doctest = true,
    clean = true,
    pages = [
        "Home" => "index.md",
        "Manual" => ["manual/saddlesystems.md",
                     "manual/timemarching.md",
                     "manual/methods.md"
                     ]
        #"Internals" => [ "internals/properties.md"]
    ],
    #format = Documenter.HTML(assets = ["assets/custom.css"])
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        mathengine = MathJax(Dict(
            :TeX => Dict(
                :equationNumbers => Dict(:autoNumber => "AMS"),
                :Macros => Dict()
            )
        ))
    ),
    #assets = ["assets/custom.css"],
    #strict = true
)


#if "DOCUMENTER_KEY" in keys(ENV)
deploydocs(
     repo = "github.com/JuliaIBPM/ConstrainedSystems.jl.git",
     target = "build",
     deps = nothing,
     make = nothing
     #versions = "v^"
)
#end
