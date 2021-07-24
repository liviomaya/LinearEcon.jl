module LinearEcon

export Model, VARModel, moving_average, var_decomp, irf, simulate, cov, corr, cholesky

using LinearAlgebra
using Distributions
using Plots
gr()
include("structures.jl")
include("modelsolution.jl")
include("timeseries.jl")

end 
