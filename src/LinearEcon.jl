module LinearEcon

export Model, Solution, model, solution, covariance, correlation, path, irf, vardecomp

using LinearAlgebra
using Plots
pyplot()
include("structures.jl")
include("auxiliary.jl")
include("treat.jl")
include("main.jl")
include("studies.jl")

end 
