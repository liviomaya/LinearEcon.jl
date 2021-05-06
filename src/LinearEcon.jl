module LinearEcon

export Model, VARModel, MovingAverage, VarDecomp, IRF, Simulation, Covariance, Correlation, CholDecomp

using LinearAlgebra
using Distributions
using Plots
gr()
include("structures.jl")
include("modelsolution.jl")
include("timeseries.jl")

end 
