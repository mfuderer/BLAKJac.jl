# Interface functions between BLAKJac optimization and underlying BLAKJac analysis functions
# Added 2022-07-11


# Since the optimizer wants to interface to a function accepting only floats, but optimization of complex values may be needed,
# this function converts a complex vector into a real vector of doubled length
function ConvertCpxToConcatenatedFloat(RFC::Vector{ComplexF64})
    length = size(RFC, 1)
    RF = zeros(length * 2)
    for i = 1:length
        RF[i] = real(RFC[i])
        RF[i+length] = imag(RFC[i])
    end
    return RF
end

# If the vector-to-be-optimized is already non-complex float, the "Unwrap" is a dummy function
function ConvertCpxToConcatenatedFloat(v::Vector{Float64})
    return v
end

# "Wrapping" is the reverse function: converting a double-length float vector into a complex vector
function ConvertConcatenatedFloatToCpx(RF::Vector{Float64})
    halfsize = size(RF, 1) ÷ 2
    RFC = zeros(ComplexF64, halfsize)
    for i = 1:halfsize
        RFC[i] = complex(RF[i], RF[i+halfsize])
    end
    return RFC
end

# Interface between (A) an optimizer that wants to optimize on a limited set FLOAT of values and (B) an analyzer that expects a full-length COMPLEX array
function ExpandCpxAndAnalyze(RF::Vector{Float64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    # convert array of floats into complex
    RFC = ConvertConcatenatedFloatToCpx(RF) # Combine float-array of double length into omplex array
    return ExpandAndAnalyze(RFC, trajectorySet, options)
end

# Interface between (A) an optimizer that wants to optimize on a limited set of values and (B) an analyzer that expects a full-length array
function ExpandAndAnalyze(RFshortPDD::Array, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    # interpolate short RF array to full length
    nTR = length(trajectorySet)
    RFdeg = ExpandRF(RFshortPDD, nTR, options)
    if ndims(RFshortPDD) == 1
        # The logic here is: if RFshort has two columns, then the second serves for a 2nd derivative of the phase;
        #   in that case it is pointless to expect complex values in the first column, which are the flip angles
        RFdeg = complex(RFdeg)  # convert to complex - if not so already
    else
        # no activity; supposed be cpx by output of ExpandRF 
    end
    penalty = BLAKJac_criterion(RFdeg, trajectorySet, options)
    return penalty
end


# Interfaces between (A) an optimizer that only wants to optimize on a FLOAT portion of a vector and (B) an analyzer that considers the full COMPLEX length
function PlugInCpxAndAnalyze(RFpart::Vector{Float64}, portionRange, RFdeg::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    RFC = ConvertConcatenatedFloatToCpx(RFpart) # Combine float-array of double length into omplex array
    return PlugInAndAnalyze(RFC, portionRange, RFdeg, trajectorySet, options)
end

# Interfaces between (A) an optimizer that only wants to optimize on a portion of a vector and (B) an analyzer that considers the full length
function PlugInAndAnalyze(RFpart::Vector{T}, portionRange, RFdeg::Vector{T}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict) where {T}
    # plug the part into the full array
    RFdeg[portionRange] = RFpart
    RFdegC = complex(RFdeg)
    penalty = BLAKJac_criterion(RFdegC, trajectorySet, options)
    return penalty
end

# Interpolates to full length.
# If there is a second column, this is interpreted as the 2nd derivative of the phase, which may be either to-be-optimized or preset
function ExpandRF(RFinPDD::Array, nTR, options)
    ph = options["opt_imposed_2nd_derivative_of_phase"]

    if (size(RFinPDD, 2) == 1)
        RFall = IncreaseResolution(RFinPDD, nTR)
        return RFall
    else
        # assumed to be one column for amplitude, other column for 2nd derivative of phase
        RFdegC = zeros(ComplexF64, nTR)
        RFallPDD = IncreaseResolution(RFinPDD, nTR)

        # if provided, overrule 'optimized' phase by predetermined (2nd derivative of) phase and convert to deciDegrees
        if (length(ph) > 0)
            RFallPDD[:, 2] .= 10.0 .* ph
        end

        # double-integrate second column and add it as a phase
        dp = 0.0
        p = 0.0
        for i in 1:nTR
            dp += deg2rad(0.1 * RFallPDD[i, 2])  # convert phase from deciDegrees to radians  
            p += dp
            RFdegC[i] = RFallPDD[i, 1] * exp(im * p)
        end
        return RFdegC
    end
end

# For a given complex vector of RF angles and a trajectory, it outputs a criterion value, according to the selected option.
# It may also, conditionally, plot an output
function BLAKJac_criterion(RFdegC::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)

    # --------------------------------------------------
    # Analyze the sequence
    noises, infocon, b1s, CSF_penalty = BLAKJac_analysis!(RFdegC, trajectorySet, options)

    # --------------------------------------------------
    # Select the criterion value ("out") given the options
    noisePenalties = Dict()
    b1Penalties = Dict()
    noisePenalties["rho"] = noises[1]
    noisePenalties["T1"] = noises[2]
    noisePenalties["T2"] = noises[3]
    noisePenalties["mean"] = mean(noises)
    noisePenalties["max"] = maximum(noises[2:3])
    noisePenalties["weighted"] = noises[2] / options["T1ref"] + noises[3] / options["T2ref"]
    b1Penalties["mean"] = mean(abs.(b1s))
    b1Penalties["max"] = maximum(abs.(b1s[2:3]))
    b1Penalties["T1"] = abs(b1s[2])
    b1Penalties["T2"] = abs(b1s[3])

    noisePenalty = noisePenalties[options["opt_focus"]]
    b1Penalty = b1Penalties[options["opt_focus"]]
    maxFlip = maximum(abs.(RFdegC))
    lowFlipFactor = 1.0 / sqrt(1.57 - 0.57 * min(maxFlip / 180.0, 2.0)) # see logbook dd. 2021-05-25  # modified 2021-12-22
    nplff = lowFlipFactor * noisePenalty
    sarLevel = mean(abs2.(RFdegC)) / options["TR"]
    sarPenalty = (sarLevel / options["sar_limit"])^10

    # --------------------------------------------------
    # Fill in default values for legacy reasons (if not provided)
    account_SAR = get(options, "account_SAR", false)
    lambda_B1 = get(options, "lambda_B1", 10.0)
    lambda_CSF = get(options, "lambda_CSF", 0.0)

    criterion = options["opt_criterion"]

    out = 0.0
    if criterion == "information content"
        out = -0.01 * infocon  # negative, since the optimizer tries to minimize the output
    # the multiplication by 0.01 serves to bring the figure to the same order of magnitude as the noise values,
    # since the NelderMead optimizer applies absolute tolerance criteria
    elseif criterion == "information content penalized"
        out = -0.01 * infocon * sqrt(1.0 / (1.3 * lowFlipFactor)) # very houtje-touwtje; the 1.3 has been chosen to make a max of approx 80 deg 'neutral' 
    else # else clause handles all situations where the criterion is based on noise, possibly combined with another criterion
        # ---------------------------------------------------------- Legacy translation 
        if criterion == "B1 sensitivity/noise mix, sar-limited"
            # legacy translation
            criterion = "noise_level"
            account_SAR = true
            lambda_B1 = 10.0
        end
        if criterion == "sar-limited noise"
            criterion = "noise_level"
            account_SAR = true
        end
        if criterion == "B1 sensitivity sar-limited"
            throw("B1 sensitivity sar-limited is not supported anymore")
        end
        # ---------------------------------------------------------- End of legacy translation

        if criterion == "noise_level"
            out = noisePenalty
        elseif criterion == "low-flip corrected noise"
            out = nplff
        else
            throw("unknown optimization criterion")
        end

        #out = out + lambda_B1*b1Penalty + lambda_CSF*CSF_penalty + account_SAR*sarPenalty
        # Modification: see justification in logbook dd. 2024-05-28
        out = sqrt(out^2 + (lambda_B1 * b1Penalty)^2 + (lambda_CSF * CSF_penalty)^2 + account_SAR * sarPenalty^2)

        negFlip = 0.0
        imFlip = 0.0
        if options["opt_keep_positive"]
            negFlip = maximum(max.(-real.(RFdegC), 0.0))
            out += negFlip
        end
    end

    # --------------------------------------------------
    # Occasionally, display intermediate output
    options["optcount"] += 1
    i = options["optcount"]
    emergeCriterion = options["opt_emergeCriterion"] # (1*1000^3) ÷ (length(RFdegC)^2)
    emerge = (i % emergeCriterion == 1)
    if emerge
        options["plotfuncs"]["close"]()
        options["plotfuncs"]["intermediate"](RFdegC, options)
        @printf("%s\n", options["stage_text"])
        @printf("Iter Noise   SAR   Sign   B1     CSF    out \n")
    end
    if (i % (emergeCriterion ÷ 5 + 1) == 1)
        # stageText = options["stage_text"]
        # @show noisePenalty, nplff, infocon, out, stageText
        # @show noisePenalty, sarPenalty, sarLevel, negFlip, out
        @printf("%4d %6.4f %6.4f %6.4f %6.4f %6.4f %8.5f\n", i, noisePenalty, sarPenalty, negFlip, b1Penalty, CSF_penalty, out)
    end

    # (Code by Hanna has been refactored into BLAKJac_analysis)

    return out
end

# Interfaces between (A) an outside world that wants to see a potentially complex array optimized using specific options
#                and (B) an optimizer that only wants to see a function accepting a float array and nothing else
function WrappedPortionOptimize(RFpart::Vector{ComplexF64}, portionRange, RFdeg::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    RFfloat = ConvertCpxToConcatenatedFloat(RFpart)
    opt_method = options["opt_method"]
    optpars = options["optpars"]

    PlugInCpxAndAnalyze_(y) = PlugInCpxAndAnalyze(y, portionRange, RFdeg, trajectorySet, options)
    res = optimize(PlugInCpxAndAnalyze_, RFfloat, opt_method(), optpars)
    RFportion = ConvertConcatenatedFloatToCpx(Optim.minimizer(res))
    return res, RFportion
end

function WrappedPortionOptimize(RFpart::Vector{Float64}, portionRange, RFdeg::Vector{Float64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    opt_method = options["opt_method"]
    optpars = options["optpars"]
    PlugInAndAnalyze_(y) = PlugInAndAnalyze(y, portionRange, RFdeg, trajectorySet, options)
    res = optimize(PlugInAndAnalyze_, RFpart, opt_method(), optpars)
    return res, Optim.minimizer(res)
end

function WrappedLowResOptimize(RFshort::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    RFfloat = ConvertCpxToConcatenatedFloat(RFshort) # Split complex array into float-array of double length
    opt_method = options["opt_method"]
    optpars = options["optpars"]
    ExpandCpxAndAnalyze_(RF) = ExpandCpxAndAnalyze(RF, trajectorySet, options)
    res = optimize(ExpandCpxAndAnalyze_, RFfloat, opt_method(), optpars)
    RFold = ConvertConcatenatedFloatToCpx(Optim.minimizer(res))
    return res, RFold
end

function WrappedLowResOptimize(RFshortPDD::Array{Float64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    opt_method = options["opt_method"]
    optpars = options["optpars"]
    ExpandAndAnalyze_(y) = ExpandAndAnalyze(y, trajectorySet, options)
    res = optimize(ExpandAndAnalyze_, RFshortPDD, opt_method(), optpars)
    return res, Optim.minimizer(res)
end



