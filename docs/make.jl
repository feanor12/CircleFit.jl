using Documenter, CircleFit

makedocs(sitename="CircleFit.jl",
	 pages = [ 
				"Overview" => "index.md"
			  "API" => "api.md"
		   "Examples" => "examples.md"
		  ])

deploydocs(
    repo = "github.com/feanor12/CircleFit.jl.git",
)
