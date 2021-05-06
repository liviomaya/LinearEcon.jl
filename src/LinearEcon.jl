module LinearEcon

export Model, VARModel, model #, solution, covariance, correlation, path, irf, vardecomp

using LinearAlgebra
using Plots
gr()
include("structures.jl")
include("modelsolution.jl")
#include("studies.jl")

end 
