# Introduced 2022-07-14
# Serves to fill in meaningful defaults for the options needed for BLAKJac optimization

function nTR(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})
    return length(trajectorySet)
end

function nky(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})
    ky = [s[i].ky for s in trajectorySet for i in 1:length(s)]
    return round(Int64, maximum(ky) - minimum(ky) + 1.0)
end
function nkz(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})
    kz = [s[i].kz for s in trajectorySet for i in 1:length(s)]
    return round(Int64, maximum(kz) - minimum(kz) + 1.0)
end

function conFill!(opt::Dict, s::String, v)
    if !(s in keys(opt))
        opt[s] = v
    end
end


"""
    BLAKJac_defaults!(trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)

Initializes the elements of the dictionary 'options' (if not already filled in) to meaningful default values.

it requires a trajectorySet as input (see TrajectorySet()), which defines the number of profiles (nTR) as well as the ky and kz ranges.
"""
function BLAKJac_defaults!(trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, options::Dict)
    conFill!(options, "rfFile", "")            # name of the file to read the RF shape from
    conFill!(options, "rflabel", "")           # For graph labels
    vnTR = nTR(trajectorySet)
    conFill!(options, "nTR", vnTR)             # number of simulated time points            
    conFill!(options, "TR", 0.01)              # in seconds     
    conFill!(options, "T1ref", 0.67)           # The reference T1 value - by now, mainly used for noise scaling and for "relative" noise deviations
    conFill!(options, "T2ref", 0.076)
    conFill!(options, "startstate", 1.0)       # Starting state of the z-magnetisation; +1 for no prepulse, -1 for inversion
    conFill!(options, "maxstate", 64)          # Length of history taken into account in simulating magnetisation from sequence
    vnky = nky(trajectorySet)
    conFill!(options, "nky", vnky)             # Number of different encoding values simulated; nTR/nky/nkz is number of samples per encoding    
    vnkz = nkz(trajectorySet)
    conFill!(options, "nkz", vnkz)
    v = round(Int64, 4 * vnTR / (vnky * vnkz))
    conFill!(options, "maxMeas", v)             # Maximum number of measurements per phase-encoding value
    conFill!(options, "handleB1", "no")         # "no", "sensitivity", "co-reconstruct"
    conFill!(options, "considerCyclic", false)   # If true, then the RF pattern (and only the RF pattern) is considered cyclic
    conFill!(options, "useSymmetry", false)     # Assume real-ness of rho, T1 and T2 and therefor symmetry in k-space
    v = [200.0^(-2), 200.0^(-2), 200.0^(-2)]
    conFill!(options, "invregval", v)           # The inverses of the regularization matrix
    conFill!(options, "useSurrogate", false)    # Using the surrogate model tends to be beneficial if at least 24 (T1,T2)-combinations are to be evaluated
    conFill!(options, "sigma_ref", 0.2)         # Reference normalized noise level for information-content estimation. See logbook around 2021-06-03
    conFill!(options, "sar_limit", 40^2 / 0.01)   # SAR limit. Initially set to an RMS level of 40 degrees at TR=10ms
    v = [(0.8183, 0.0509), (0.67, 0.2066), (1.2208, 0.1253), (0.2465, 0.0461), (2.2245, 0.3082), (0.3677, 0.0461), (0.3677, 0.1253)]
    conFill!(options, "T1T2set", v)        # set of (T1,T2) combinations (in the log(T/Tref) space) for which n-factor and a-factor will be evaluated
    conFill!(options, "sizeSteps", [10])        # In the first stage of the optimization procedure, the resolutions of the sequence that will be optimized upon

    conFill!(options, "portion_size", 50)       # In an optimization procedure, the size of the sequence-portion that will be optimized "in-context" 
    conFill!(options, "portion_step", 30)       # The portion-to-be-optimized is stepped over this number of TRs
    v = length(options["sizeSteps"])
    conFill!(options, "opt_iterations_limit", v) # Possibility to cut short the optimization process, e.g. to first stage only 
    v = Optim.Options(time_limit=3000.0, iterations=100000, f_tol=1.0e-2, g_tol=1.0e-3)
    conFill!(options, "optpars", v)             # Parameters for the NelderMead optimization algorithm
    conFill!(options, "opt_method", NelderMead)
    conFill!(options, "opt_expand_type", "spline") # in ["piecewise", "spline", "nodes"]
    conFill!(options, "opt_keep_positive", false) # if true, then disallow negative RF values
    conFill!(options, "opt_initialize", "cRandom30") # type of initialization for optimizer, 
    # in ["cRandom30", "ones", "ernst", "init_angle", "quadraticPhase30", "RF_shape"]
    conFill!(options, "opt_complex", true)      # complex optimization, each excitation having its own phase
    conFill!(options, "opt_slow_phase", false)  # allow quadratic phase increment
    conFill!(options, "opt_imposed_2nd_derivative_of_phase", []) # allows to provide a preselected 2nd derivative op the phase, in radiants
    conFill!(options, "opt_focus", "mean")      # to focus on one of ["rho","T1","T2","mean","weighted","max"]
    conFill!(options, "emphasize_low_freq", false) # if switched on, then the optimization will steer more on reducing the noise factor on low spatial frequencies
    conFill!(options, "opt_criterion", "sar-limited noise") # "noise_level", "low-flip corrected noise", "information content", "information content penalized"
    conFill!(options, "account_SAR", true)      # if true, then the SAR is taken into account
    conFill!(options, "B1metric", "multi_point")# "multi_point" for optimization around B1=1.0; "derivative_at_1" for analysis
    conFill!(options, "lambda_B1", 10.0)        # penalty for B1 sensitivity
    conFill!(options, "lambda_CSF", 0.0)        # penalty for CSF sensitivity
    conFill!(options, "opt_account_maxFlip", true) # take into account that RF pulses need time
    conFill!(options, "opt_emergeCriterion", 1000) # how seldomly do we want an in-between result?
    conFill!(options, "plot", false)             # is plotting of intermediate results required?
    conFill!(options, "plottypes", [])           # collection of "first", "trajectories", "weights", "angles", "noisebars", "noiseSpectrum", "infocon"
    conFill!(options, "optcount", 0)
    conFill!(options, "stage_text", "")

    plotfuncs = setPlotFuncs()
    conFill!(options, "plotfuncs", plotfuncs)
end

