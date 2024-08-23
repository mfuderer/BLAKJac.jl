# Introduced 2022-07-14
# Serves to fill in meaningful defaults for the options needed for BLAKJac optimization
"""
    BLAKJac_defaults!(trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, options::Dict)

Initializes the elements of the dictionary 'options' (if not already filled in) to meaningful default values. 
    
A trajectorySet is required as input (see the `TrajectorySet` function), which defines the number of profiles (nTR) as well as the ky and kz ranges.
"""
function BLAKJac_defaults!(trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, options::Dict)
    get!(options, "rfFile", "")            # name of the file to read the RF shape from
    get!(options, "rflabel", "")           # For graph labels
    vnTR = nTR(trajectorySet)
    get!(options, "nTR", vnTR)             # number of simulated time points            
    get!(options, "TR", 0.01)              # in seconds     
    get!(options, "T1ref", 0.67)           # The reference T1 value - by now, mainly used for noise scaling and for "relative" noise deviations
    get!(options, "T2ref", 0.076)
    get!(options, "startstate", 1.0)       # Starting state of the z-magnetisation; +1 for no prepulse, -1 for inversion
    get!(options, "maxstate", 64)          # Length of history taken into account in simulating magnetisation from sequence
    vnky = nky(trajectorySet)
    get!(options, "nky", vnky)             # Number of different encoding values simulated; nTR/nky/nkz is number of samples per encoding    
    vnkz = nkz(trajectorySet)
    get!(options, "nkz", vnkz)
    v = round(Int64, 4 * vnTR / (vnky * vnkz))
    get!(options, "maxMeas", v)             # Maximum number of measurements per phase-encoding value
    get!(options, "handleB1", "no")         # "no", "sensitivity", "co-reconstruct"
    get!(options, "considerCyclic", false)   # If true, then the RF pattern (and only the RF pattern) is considered cyclic
    get!(options, "useSymmetry", false)     # Assume real-ness of rho, T1 and T2 and therefor symmetry in k-space
    v = [200.0^(-2), 200.0^(-2), 200.0^(-2)]
    get!(options, "invregval", v)           # The inverses of the regularization matrix
    get!(options, "useSurrogate", false)    # Using the surrogate model tends to be beneficial if at least 24 (T1,T2)-combinations are to be evaluated
    get!(options, "sigma_ref", 0.2)         # Reference normalized noise level for information-content estimation. See logbook around 2021-06-03
    get!(options, "sar_limit", 40^2 / 0.01)   # SAR limit. Initially set to an RMS level of 40 degrees at TR=10ms
    v = [(0.8183, 0.0509), (0.67, 0.2066), (1.2208, 0.1253), (0.2465, 0.0461), (2.2245, 0.3082), (0.3677, 0.0461), (0.3677, 0.1253)]
    get!(options, "T1T2set", v)        # set of (T1,T2) combinations (in the log(T/Tref) space) for which n-factor and a-factor will be evaluated
    get!(options, "sizeSteps", [10])        # In the first stage of the optimization procedure, the resolutions of the sequence that will be optimized upon

    get!(options, "portion_size", 50)       # In an optimization procedure, the size of the sequence-portion that will be optimized "in-context" 
    get!(options, "portion_step", 30)       # The portion-to-be-optimized is stepped over this number of TRs
    v = length(options["sizeSteps"])
    get!(options, "opt_iterations_limit", v) # Possibility to cut short the optimization process, e.g. to first stage only 
    v = Optim.Options(time_limit=3000.0, iterations=100000, f_tol=1.0e-2, g_tol=1.0e-3)
    get!(options, "optpars", v)             # Parameters for the NelderMead optimization algorithm
    get!(options, "opt_method", NelderMead)
    get!(options, "opt_expand_type", "spline") # in ["piecewise", "spline", "nodes"]
    get!(options, "opt_keep_positive", false) # if true, then disallow negative RF values
    get!(options, "opt_initialize", "cRandom30") # type of initialization for optimizer, 
    # in ["cRandom30", "ones", "ernst", "init_angle", "quadraticPhase30", "RF_shape"]
    get!(options, "opt_complex", true)      # complex optimization, each excitation having its own phase
    get!(options, "opt_slow_phase", false)  # allow quadratic phase increment
    get!(options, "opt_imposed_2nd_derivative_of_phase", []) # allows to provide a preselected 2nd derivative op the phase, in radiants
    get!(options, "opt_focus", "mean")      # to focus on one of ["rho","T1","T2","mean","weighted","max"]
    get!(options, "emphasize_low_freq", false) # if switched on, then the optimization will steer more on reducing the noise factor on low spatial frequencies
    get!(options, "opt_criterion", "sar-limited noise") # "noise_level", "low-flip corrected noise", "information content", "information content penalized"
    get!(options, "account_SAR", true)      # if true, then the SAR is taken into account
    get!(options, "B1metric", "multi_point")# "multi_point" for optimization around B1=1.0; "derivative_at_1" for analysis
    get!(options, "lambda_B1", 10.0)        # penalty for B1 sensitivity
    get!(options, "lambda_CSF", 0.0)        # penalty for CSF sensitivity
    get!(options, "opt_account_maxFlip", true) # take into account that RF pulses need time
    get!(options, "opt_emergeCriterion", 1000) # how seldomly do we want an in-between result?
    get!(options, "plot", false)             # is plotting of intermediate results required?
    get!(options, "plottypes", [])           # collection of "first", "trajectories", "weights", "angles", "noisebars", "noiseSpectrum", "infocon"
    get!(options, "optcount", 0)
    get!(options, "stage_text", "")

    plotfuncs = setPlotFuncs()
    get!(options, "plotfuncs", plotfuncs)
end

