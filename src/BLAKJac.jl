module BLAKJac

# Load external packages
using BlochSimulators
using CUDA
using DelimitedFiles
using Distributed
using DistributedArrays
using FileIO
using Interpolations
using LinearAlgebra
using Optim
using Printf
using Random
using StaticArrays
using Statistics
using StructArrays

include("analysis_structs.jl")
include("utils.jl")
include("analysis.jl")
include("defaults.jl")
include("interfaceOptToAn.jl")
include("optimize.jl")

# Plotting functions: If PyPlot is available in the environment from which BLAKJac is loaded, then PyPlot-based plotting functions are loaded. Otherwise nothing is plotted.
# Users can define their own plotting functions and load them in the same way as PyPlot is loaded.
include("../plot/set_plotfuncs.jl")
include("../plot/pyplot.jl")

export BLAKJac_analysis!, BLAKJac_criterion, WrappedLowResOptimize, WrappedPortionOptimize
export BLAKJac_optimize, BLAKJac_defaults!, TrajectorySet, PlotAmplAndPhaseDeriv

end