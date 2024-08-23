
# BLAKJac_optimize: optimization procedures according to BLAKJac
#  added 2022-07-13 to MRSTAT

"""
    BLAKJac_optimize(trajectorySet, options::Dict, mersenneTwisterSeed=1)

Optimization procedure according to BLAKJac

For a given ky pattern (or (ky,kz)-pattern) and given some options, it generates a sequence of RF pulses pretending to be optimized.
The seed parameter can be provided as an integer, to allow different realizations of random-value initializations.
(By default, it should reproducably generate the same sequence on every call)

A meaningful set of options can be pre-defined by "options = Dict(); BLAKJac.BLAKJac_defaults!(trajectorySet, options)"

To understand some option parameters, one must know that the optimization consists of two stages (of which often only one is used);
the first stage optimizes cubic-spline-interpolated smooth RF-patterns, with ever increasing number of degrees of freedom ("increasing resolution").
In the second stage, successive portions of the sequence are "fine-tuned" by optimizing the sequence portion by portion, taking all of the optimized 
"past" for granted

Option parameters are:
  - `nky`             Extent of ky
  - `sizeSteps`       Array of integers; provides subsequent "resolutions" in the first stage
  - `portion_size`    The size of a portion to be optimized
  - `portion_step`    after one portion is considered sufficiently optimized, the "next" portion is addressed, which is `portion_step` further on.
                      Preferably, `portion_step` is smaller than `portion_size`.
  - `opt_iterations_limit`  Int. Controls the total of first-stage and second-stage steps. To switch off the second stage, it has to be length of sizeSteps.
                        Setting it to an arbitrary high value means "full processing"
  - `opt_complex`     Bool. Allow all pulses to have a phase. Deprecated.
  - `opt_slow_phase`  Bool. A very controlled way to optimize the phase: what is optimized is a low-resolution of the second derivative of the phase
  - `opt_imposed_2nd_derivative_of_phase`        (Vector over nTR elements) To be filled in if the 2nd derivative of phase is not to be optimized - 
                                              i.e. if it is pre-defined
  - `opt_initialize`  Controls how the RF-pattern is optimized before entering optimization. Can be one of `"ones"` (all angles are 1 degree), 
                        `"ernst"` (all angles are set to the Ernst angle given a reference T1 and the TR), `"init_angle"` (all set to `opt_init_angle`),
                        `"cRandom30"` (put it to random values with a deviation of 30 degrees) or `"RF_shape"` (initial shape provided by option 
                        `"rfFunction"`)
  - `rfFunction`      (see above)
  - `opt_init_angle`  (see above)
  - `opt_focus`       Which property should the optimizer focus on? Can be any of `"rho"`, `"T1"`, `"T2"`, `"mean"`, `"max"` or `"weighted"`; 
                        the last one 
                      combines the *relative* noise levels of T1 and T2 (relative to T1ref resp. T2ref)
  - `opt_criterion`   Can be any of `"noise_level"`, `"low-flip corrected noise"`, `"sar-limited noise"`, `"information content"`,
                     `"information content penalized"`, 
                      `"B1 sensitivity sar-limited"`, `"B1 sensitivity/noise mix, sar-limited"`
  - `lambda_B1`       A regularization parameter for the B1 sensitivity
  - `lambda_CSF`      A regularization parameter for the CSF sensitivity
  - `opt_keep_positive`   Bool. If true, negative angles are heavily penalized
  - `opt_emergeCriterion` Int. Calls for intermediate displays of the RF shape every opt_emergeCriterion iterations
  - `opt_method`      optimizer function, e.g. NelderMead
  - `optpars`         optimization parameters - see Optim package.
"""
function BLAKJac_optimize(trajectorySet, options::Dict, mersenneTwisterSeed=1)
    rSeed = MersenneTwister(mersenneTwisterSeed)

    nTR = length(trajectorySet)
    portionSize = options["portion_size"]
    portionStep = options["portion_step"]
    options["plot"] = false                     # No plotting during iterations
    options["plottypes"] = []                      # No plotting during iterations
    cutShort = options["opt_iterations_limit"]
    sizeSteps = options["sizeSteps"]

    startTime = time()

    firstPhase = length(sizeSteps)
    if (options["opt_expand_type"] == "spline")

        RFoldPDD = InitializeRFPattern(rSeed, nTR, options)  # "PDD" = "Phase in DeciDegrees"
        # If we want to optimize a 2nd derivative of the phase, this will be stored in the 2nd column of RFoldPDD,
        # and the units thereof will be in deciDegrees (actually deciDegrees per TR^2). 
        # The choice for such an apparently weird unit is that the two columns have similar ranges, which
        # is beneficial for the optimization process.

        # first stage: gradually increasing resolution
        for iterCount = 1:min(length(sizeSteps), cutShort)
            thisSize = sizeSteps[iterCount]
            options["stage_text"] = "first phase, size is $thisSize"
            RFshortPDD = IncreaseResolution(RFoldPDD, thisSize)

            res, RFoldPDD = BLAKJac.WrappedLowResOptimize(RFshortPDD, trajectorySet, options)

            @show time() - startTime
            @show res
            # Main.PyPlot.pause(2.0)
            sleep(2.0)
        end
        # interpolate short RF array to full length
        RFdeg = ExpandRF(RFoldPDD, nTR, options)
    else
        @assert(false, "unsuported expand type")
    end

    # Second stage: optimize on successive portions of the sequence
    portionStart = 1
    iterCount = firstPhase + 1
    while ((portionStart <= nTR) && (iterCount <= cutShort) && (options["opt_slow_phase"] == false))
        portionEnd = min(nTR, portionStart + portionSize - 1)
        portionRange = portionStart:portionEnd
        options["stage_text"] = "second phase, considered range is $portionRange"

        RFpart = RFdeg[portionRange]
        res, RFportion = BLAKJac.WrappedPortionOptimize(RFpart, portionRange, RFdeg, trajectorySet, options)
        RFdeg[portionRange] = RFportion

        FoldToBaseRange!(RFdeg)

        @show time() - startTime
        @show res
        pause(2.0)

        portionStart += portionStep
        iterCount += 1
    end

    return RFdeg # ,RFold
end

# Generate start pattern for the optimization process according to the options
#  (This overload exists for legacy reasons and can be called if the 2nd-derivative-of-phase is not involved)
function InitializeRFPattern(rSeed::MersenneTwister, nTR, options, returnCpx::Bool)
    sizeSteps = options["sizeSteps"]
    initSize = sizeSteps[1]
    RFold = ones(ComplexF64, initSize)
    type = options["opt_initialize"]
    if type == "ones"
        # already filled in
    elseif type == "ernst"
        ernst = acosd(exp(-options["TR"] / options["T1ref"]))
        RFold .*= ernst
    elseif type == "init_angle"
        RFold .*= options["opt_init_angle"]
    elseif type == "cRandom30"
        RFold = 30.0 * randn(rSeed, ComplexF64, initSize)
    elseif type == "quadraticPhase30"
        half = (initSize - 1.0) / 2.0
        RFold = [30.0 * exp(im * pi * ((i - 1 - half) / half)^2) for i in 1:initSize]
    elseif type == "RF_shape"
        rfFunction = options["rfFunction"]
        RFpredef = rfFunction(nTR, options["nky"])
        stride = nTR / sizeSteps[1]
        RFold = [RFpredef[convert(Int64, round((i - 0.5) * stride + 0.5))] for i in 1:sizeSteps[1]]
    else
        @assert(false)
    end

    if returnCpx
        return RFold
    else
        return real.(RFold)
    end
end

# Generate start pattern for the optimization process according to the options
#    (generic overload)
function InitializeRFPattern(rSeed::MersenneTwister, nTR, options::Dict)
    if options["opt_complex"]
        return InitializeRFPattern(rSeed, nTR, options, true)  # returns 1D array of complex values
    else
        if options["opt_slow_phase"]
            RFdeg = InitializeRFPattern(rSeed, nTR, options, false)
            RFall = zeros(length(RFdeg), 2)

            type = options["opt_initialize"]
            if type == "ones"
                RFall = ones(length(RFdeg), 2)
            elseif type == "ernst"
                # keep the phase'' to zero
            elseif type == "init_angle"
                RFall .+= options["opt_init_phase_angle"]
            elseif type == "cRandom30"
                randn!(rSeed, RFall) # random values around 0, with 1 degree deviation
            else
                error("Initialization type not supported with small quadratic phase")
            end

            RFallPDD = 10.0 .* RFall        # convert phase in deciDegrees ("PDD" = "Phase in DeciDegrees")
            RFallPDD[:, 1] = RFdeg           # first column to contain amplitude (deg), 
            #   2nd column to contain 2nd derivative of phase (deciDegrees)
            return RFallPDD                 # returns 2D array of floats
        else
            return InitializeRFPattern(rSeed, nTR, options, false) # returns 1D array of floats
        end
    end
end

# Bring elements to(wards) the (-180,+180) phase range
function FoldToBaseRange!(RF)
    for (i, v) in enumerate(RF)
        if (abs(v) > 180.0)
            wrapvalue = 360.0 * v / abs(v)
            RF[i] -= wrapvalue
        end
    end
end

function IncreaseResolution(RFold, thisSize)
    RFshort = similar(RFold, thisSize, size(RFold, 2))
    if size(RFold, 1) > 1
        # interpolate
        for c = 1:size(RFold, 2)
            itp = Interpolations.interpolate(RFold[:, c], BSpline(Cubic(Natural(OnGrid())))) # InPlace,Natural
            RFshort[:, c] = itp.(LinRange(1, size(RFold, 1), thisSize))          # increase resolution to next level 
        end
    else
        # degenerate case: copy single-element
        for c = 1:size(RFold, 2)
            RFshort[:, c] .= RFold[1, c]
        end
    end
    if size(RFold, 2) == 1                    # somewhat clumsy way of removing singleton
        RFshort = RFshort[:, 1]
    end
    return RFshort
end