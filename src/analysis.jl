
"""
BLAKJac_analysis!(RFdeg::Vector{ComplexF64}, trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, options::Dict, saved_H::Dict=Dict())

BLAKJac_Analysis! predicts noise levels, given an RF pattern and a phase encoding pattern.

# Inputs: 
RF pattern,   phase-encoding pattern, and a bunch of parameters (see below)

# Output:
- noisesAll   An array of 3 noise values (rho, T1, T2)  
- ItotAll     Information content (still a very raw metric)
- b1factorsRMS An array of couplings between (rho, T1, T2) and B1 (and that RMS over the T1T2set); only calculated if handleB1=="sensitivity"

# InOut:
A dictionary of H matrix arrays (labeled by combination of T1T2-probing index and RF pattern name)

# Parameters:
(Note: A meaningful set of options can be pre-defined by "options = Dict(); BLAKJac.BLAKJac_defaults!(trajectorySet, options)")

- `account_SAR`   If true, the SAR is taken into account in the optimization
- `considerCyclic`   If set to true, it is assumed that the magnetization state at the beginning of the sequence is the same as at the end thereof.
- `emphasize_low_freq`   If true, a higher weight is given to the most central region of (ky,kz)
- `handleB1`   By default ("no"), it is assumed that the only parameters to be reconstructed are T1, T2 and rho; if handleB1=="co-reconstruct",
            then BLAKJac assumes that this is a reconstructable parameter as well. If handleB1=="sensitivity", then calculate b1factorRMS
- `invregval`   An array for (inverses of) the diagonals of a regularization matrix assumed in reconstruction:
            A value of 0 assumes no imposed regularisation; a value of 1 assumes that the reconstruction very strongly imposes that the result
            should equal the reference values for T1 and T2. The breakeven occurs when invregval is on the order of `abs_sensitivity^2`,
            where `abs_sensitivity` is on the order of 0.05, depending on T1, T2 and sequence (very roughly, T2/T1). If the expected SNR is
            much larger than 1 (e.g. 10), the invregval should be significantly smaller, e.g. 0.005^2.
- `lambda_B1`  A regularization parameter for the B1 sensitivity
- `lambda_CSF` A regularization parameter for the CSF sensitivity
- `maxMeas`   By default (BLAKJac_defaults), it is set to 4 times the maximum number of times that a specific (ky,kz)-combination will be seen during the measurement. 
- `maxstate`   the maximum number of state that the EPG-estimation of signal and derivatives will take into account
- `nky`   By default (BLAKJac_defaults), the range of applied ky values 
- `nkz`   By default (BLAKJac_defaults), the range of applied kz values 
- `plotfuncs`   (Filled in by BLAKJac_defaults!(); do not adapt)
- `plottypes`   An array of strings. If the array contains a recognized string, a plot is produced
| String | description |
| :------ | :--------- |
| "first" | plots the RF shape alongside the ky and the kz of the first sample of each trajectory |
| "trajectories" | plots the first 10 trajectories |
| "noisebars" |  a barplot of the noise levels |
| "weights" |  colorful plot of T1- and T2-sensitivities over time |
| "original_jacobian" | (to be described) |
| "noiseSpectrum" | Graphs of noise spectral density over ky, one for T1 map and one for T2 map |
| "infocon" |      (to be described) |

- `rfFile`   Text used to tag the H-matrix dictionary 
- `sigma_ref`   (used in calculating information content)
- `startstate`   either 1 (meaning no inversion prepulse) or -1 (meaning an inversion prepulse one TR prior to actual sequence)
- `T1ref`   The "reference" T1 value 
- `T2ref`   The "reference" T2 value
- `T1T2set`   a set of (T1,T2) values around which BLAKJac will be evaluated. Typically "the 7-points mix" (ToDo: redesign this)
- `TR`  in seconds
- `useSurrogate`  should be false by default; if true, it uses a polynomial estimate of signal and derivatives, rather than the actual EPG estimate
- `useSymmetry`   if true, it is assumed that the phase of proton density is very smooth, such that, for BLAKJac analysis, all output can be considerd to be real
"""
function BLAKJac_analysis!(
    RFdeg::Vector{<:Complex}, 
    trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, 
    options::Dict, 
    saved_H::Dict=Dict()
    )

    # Extract the test values for T1 and T2 from the options dictionary
    T1T2set = options["T1T2set"]
    
    # Check whether B1 is to be considered as a nuisance parameter
    nNuisances = _calculate_num_nuisance_parameters(options)

    # Assemble SPGR sequence simulator (FISP3D) 
    sequence = _assemble_FISP3D(options, RFdeg)

    # Plot RF train and trajectory (if supported plotting backend is available)
    _plot_sequence(options, RFdeg, trajectorySet)

    # Estimate the noise levels, information content and B1 sensitivity factor
    B1_coupling_factors, noise_values, information_content = BLAKJacOnT1T2set(T1T2set, sequence, trajectorySet, options, saved_H, nNuisances)

    # Calculate the CSF penalty, if any
    CSF_penalty = _calculate_csf_penalty(options, sequence, T1T2set)

    # Plot noise levels (if supported plotting backend is available)
    _plot_noise_levels(options, RFdeg, trajectorySet, noise_values)

    return noise_values, information_content, B1_coupling_factors, CSF_penalty
end

"""
In previous versions of BLAKJac, a `resource` argument was required. In the current version it is removed because it was not used in the function. For backwards-compatability purposes, we allow BLAKJac_analysis! to be called with the `resource` argument. 
"""
function BLAKJac_analysis!(resource,
    RFdeg::Vector{<:Complex}, 
    trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, 
    options::Dict, 
    saved_H::Dict=Dict()
    )

    BLAKJac_analysis!(RFdeg, trajectorySet, options, saved_H)
end

"""
    BLAKJacOnT1T2set(T1T2set, sequence, trajectorySet, options, saved_H, nNuisances)

Predicts noise levels, information content and B1 coupling factors for a given combination of RF pattern (that is contained with `sequence`) and phase encoding pattern (`trajectorySet`). 

Whereas `BLAKJacOnSingleT1T2` uses a single pair of T₁ and T₂ values, `BLAKJacOnT1T2set` bases the predictions on a set of T₁ and T₂ values.
"""
function BLAKJacOnT1T2set(T1T2set, sequence, trajectorySet, options, saved_H, nNuisances)

    # Initialize accumulators to calculate averaged result over the T1T2set
    nPars = 3 # rho, T1, T2
    b1factors2All = zeros(nPars)
    noisesAll = zeros(nPars)
    b1factorsLast = zeros(nPars)
    ItotAll = 0.0

    # Loop over all probe-values of (T1,T2)
    for (index, (T1test, T2test)) in enumerate(T1T2set)
        B1test = 1.0

        # Calculate noise levels, information content and B1 sensitivity factor
        H, noises, Itot, b1factors = BLAKJacOnSingleT1T2(T1test, T2test, B1test, nNuisances, sequence, trajectorySet, options)

        # Accumulate results
        b1factors2All .+= (b1factors) .^ 2
        b1factorsLast = b1factors
        noisesAll .+= noises
        ItotAll += Itot

        if (@isdefined saved_H)
            saved_H["$index for $(options["rfFile"])"] = H
        end
    end #for loop over probe-set of (T1,T2)

    b1factorsRms = length(T1T2set) == 1 ? b1factorsLast : sqrt.(b1factors2All ./ size(T1T2set))
    noisesAll ./= size(T1T2set)
    ItotAll /= only(size(T1T2set))

    return b1factorsRms, noisesAll, ItotAll
end

"""
    BLAKJacOnSingleT1T2(T1test, T2test, B1test, nNuisances, spgr::BlochSimulators.FISP3D, trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, options::Dict)

Predicts noise levels, information content and B1 coupling factors for a given combination of RF pattern (that is contained with `sequence`) and phase encoding pattern (`trajectorySet`) using a **single pair** of T₁ and T₂ values. 
"""
function BLAKJacOnSingleT1T2(T1test, T2test, B1test, nNuisances, spgr::BlochSimulators.FISP3D, trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, options::Dict)

    # Load options and calculated derived values
    nTR = length(trajectorySet)
    nky = options["nky"]
    nkz = options["nkz"]
    useSym = options["useSymmetry"]
    nkyEff = useSym ? nky ÷ 2 : nky
    nkyEff = max(1, nkyEff)
    nPars = 3

    note = @sprintf("for (T1,T2)=(%6.1f,%6.1f)", 1000 * T1test, 1000 * T2test)
    
    takeB1asVariable = (options["handleB1"] == "co-reconstruct")

    nParsX = nPars
    nNuisancesX = nNuisances
    if (nPars == 3) && (nNuisances == 1) && takeB1asVariable
        nParsX = 4
        nNuisancesX = 0
    end

    # Calculate local weights (magnetization and derivatives)
    wlocal = _calculate_local_weights(spgr, options, nPars, nTR, nNuisances, T1test, T2test, B1test)

    # Plot calculated weights (if supported plotting backend is available)
    _plot_weights(options, trajectorySet, wlocal)

    # Plot original Jacobian (if supported plotting backend is available)
    # TODO: They seem to plot the same thing (wlocal)?
    _plot_original_jacobian(options, wlocal)

    # analyze Jacobian
    wmat = _analyze_jacobian(nky, nkz, nkyEff, useSym, nPars, nNuisances, nTR, wlocal, trajectorySet, options)

    # analyze matrices for all ky
    sumH, H, Hdiag = _calculate_H_matrices(wmat, takeB1asVariable, nPars, nky, nkz, useSym, nParsX, nkyEff, options)

    # B1 sensitivity factor
    b1factors = _calculate_b1_sensitivity(nPars, nkz, nNuisances, nNuisancesX, options)

    # A scaling factor is introduced that should normalize the noise level to 1 in case of Rho-only reconstruction
    # Calculate "Normalized Expected Signal squared" (nes2)
    nes2 = _calculate_nes2(nTR, options)
    
    # output diagonal elements
    _plot_noise_spectrum(H, nes2, note, options)

    # output RMS of Hessian
    noises = _calculate_noises(sumH, nkyEff, nkz, nes2)
    
    # info content
    Itot = _calculate_information_content(Hdiag, options, nes2, note, nParsX)

    return H, noises, Itot, b1factors
end

############################################################################################
######################## HELPER FUNCTIONS ##################################################
############################################################################################

function _calculate_local_weights(spgr, options, nPars, nTR, nNuisances, T1test, T2test, B1test)

    parameters = (nNuisances > 0) ? T₁T₂B₁(T1test, T2test, B1test) : T₁T₂(T1test, T2test)
    parameters = StructVector([parameters])
    fit_parameters = (nNuisances > 0) ? (:T₁, :T₂, :B₁) : (:T₁, :T₂)

    wlocal = zeros(ComplexF64, nTR, nPars + nNuisances)

    useSurrogate = options["useSurrogate"]

    if (useSurrogate)
        error("Surrogate model not supported at the moment")
    end

    m = simulate_magnetization(spgr, parameters)
    ∂m = simulate_derivatives_finite_difference(fit_parameters, m, spgr, parameters)

    wlocal[:, 1] = m
    wlocal[:, 2] = ∂m.T₁ .* T1test
    wlocal[:, 3] = ∂m.T₂ .* T2test
    if nNuisances > 0
        wlocal[:, 4] = ∂m.B₁
    end
    
    return wlocal
end

function _calculate_H_matrices(wmat, takeB1asVariable, nPars, nky, nkz, useSym, nParsX, nkyEff, options)

    invRegB1 = takeB1asVariable ? 0.0 : sqrt(prevfloat(Inf)) # corretion 2024-06-04. Logic: if B1 is to be co-reconstructed, 
    # then it is not regularized; if it is to be sensitivity-analyzed, 
    # it should not influence the estimate of T1T2rho noise, so it is
    # very heavily ('infinitely') regularized
    invReg = Diagonal([options["invregval"]; [invRegB1]])

    # analyze matrices for all ky
    sumH = zeros(ComplexF64, nParsX, nParsX)
    H = zeros(ComplexF64, nkyEff, nkz, nParsX, nParsX)
    Hdiag = zeros(ComplexF64, nkyEff, nkz, nParsX)
    
    for thiskz = 1:nkz
        for thisky = 1:nkyEff
            W = wmat[thisky, thiskz, :, 1:nParsX]
            JhJ = adjoint(W) * W
            Hlocal = inv(JhJ + invReg[1:nParsX, 1:nParsX])
            H[thisky, thiskz, :, :] = Hlocal
        end
    end

    # Addition 2021-10-06, taking into account that the noise-factors on low k-values also multiply a variety of imperfections
    H_emphasized = copy(H)
    if (options["emphasize_low_freq"])
        for ky in 1:nkyEff
            for kz in 1:nkz
                kyVal = useSym ? ky - 1 : ky - nky ÷ 2 - 1.0
                kzVal = kz - nkz ÷ 2 - 1.0
                #a = 3.0; b = 20.0  # 'a' is something like 'FOV/objectsize'; b is relative importance of artefacts at k=0
                a = 3.0
                b = 3.0  # 'a' is something like 'FOV/objectsize'; b is relative importance of artefacts at k=0
                emphasis = 1 + b / ((kyVal / a)^2 + (kzVal / a)^2 + 1.0^2)
                H_emphasized[ky, kz, :, :] = H[ky, kz, :, :] .* emphasis
            end
        end
    end

    sumH = sum(H_emphasized, dims=1)
    sumH = sum(sumH, dims=2)

    for i in 1:nPars
        Hdiag[:, :, i] = H[:, :, i, i]
    end
    
    return sumH, H, Hdiag

end

function _calculate_noises(sumH, nkyEff, nkz, nes2)
    sigmaRho = sqrt(abs(sumH[1, 1, 1, 1]) / nkyEff / nkz)
    sigmaT1 = sqrt(abs(sumH[1, 1, 2, 2]) / nkyEff / nkz)
    sigmaT2 = sqrt(abs(sumH[1, 1, 3, 3]) / nkyEff / nkz)
    noises = [sigmaRho, sigmaT1, sigmaT2]
    noises .*= sqrt(nes2)
    return noises
end

function _calculate_sumH(H_emphasized)

    sumH = sum(H_emphasized, dims=1)
    sumH = sum(sumH, dims=2)

    return sumH
end

function _calculate_nes2(nTR, options)

    # nes2 = "Normalized Expected Signal squared"

    # A scaling factor is introduced that should normalize the noise level to 1 in case of Rho-only reconstruction
    # questionable whether ok if useSym
    # The extra factor of 2.0, is yet to be explained
    TR = options["TR"]
    T1ref = options["T1ref"]
    T2ref = options["T2ref"]
    nky = options["nky"]
    nkz = options["nkz"]
    nes2 = 2.0 * 4.0 * (nTR * TR + T1ref) * T2ref / (nky * nkz * (T1ref + T2ref)^2)

    return nes2
end

function _calculate_information_content(Hdiag, options, nes2, note, nParsX)

    σref2 = options["sigma_ref"]^2
    nky = options["nky"]
    nkz = options["nkz"]
    useSym = options["useSymmetry"]
    nkyEff = useSym ? nky ÷ 2 : nky
    nkyEff = max(1, nkyEff)

    kym = useSym ? 1 : nky ÷ 2 + 1
    kzm = nkz ÷ 2 + 1
    pevaly = [((i - kym)^2 + 1^2) / ((nky ÷ 2)^2 + 1^2) for i in 1:nkyEff]
    pevalz = [((i - kzm)^2 + 1^2) / ((nkz ÷ 2)^2 + 1^2) for i in 1:nkz]

    I = zeros(Float64, nkyEff, nkz, nParsX)

    if (nkz > 1)
        for k in CartesianIndices((1:nkyEff, 1:nkz))
            ps = (pevaly[k[1]] + pevalz[k[2]])^(-1.5)                  # read "signal power model"
            pn = nes2 * σref2 * abs.(Hdiag[k[1], k[2], :])
            I[k[1], k[2], :] = log.(ps ./ pn .+ 1.0)
        end
    else
        for k in 1:nkyEff
            ps = (pevaly[k])^(-1.0)                  # read "signal power model"
            pn = nes2 * σref2 * abs.(Hdiag[k, 1, :])
            I[k, 1, :] = log.(ps ./ pn .+ 1.0)
        end
    end

    if "infocon" ∈ options["plottypes"]
        options["plotfuncs"]["infocon"](I, note)
    end

    Imetric = Dict()
    Imetric["rho"] = sum(I[:, :, 1])
    Imetric["T1"] = sum(I[:, :, 2])
    Imetric["T2"] = sum(I[:, :, 3])
    Imetric["mean"] = sum(I)
    Imetric["max"] = sum(I)
    Imetric["weighted"] = sum(I)
    Itot = Imetric[options["opt_focus"]]

    return Itot
end


function _calculate_b1_sensitivity(nPars, nkz, nNuisances, nNuisancesX, options)

    b1factors = zeros(nPars)
    b1factors2 = zeros(nPars)

    if nNuisancesX > 0
        B1metric = get(options, "B1metric", "multi_point")
        if B1metric == "derivative_at_1"
            # calculate B1 sensitivity factor from W matrix
            ky0 = useSym ? 1 : nky ÷ 2 + 1
            kz0 = nkz ÷ 2 + 1
            W = wmat[ky0, kz0, :, 1:nPars]
            Wb1 = wmat[ky0, kz0, :, nPars+1]
            JhJ = adjoint(W) * W
            H0 = inv(JhJ + invReg[1:nPars, 1:nPars])
            b1fcpx = H0 * adjoint(W) * Wb1
            b1factors = abs.(b1fcpx)

            ky0 = useSym ? 1 : nky ÷ 2 + 1
            kz0 = nkz ÷ 2 + 1
            W = wmat[ky0, kz0, :, 1:nPars]
            Wb1 = wmat[ky0, kz0, :, nPars+1]

            if useSym
                # Repair hack: remove again the extra conjugate rows of the Jacobian 
                W = wmat[ky0, kz0, 1:2:end, 1:nPars]
                Wb1 = wmat[ky0, kz0, 1:2:end, nPars+1]
            end
            JhJ = adjoint(W) * W
            H0 = inv(JhJ + invReg[1:nPars, 1:nPars])

            b1fcpx = H0 * adjoint(W) * Wb1
            b1factors = real.(b1fcpx)
        elseif B1metric in ["multi_point", "multi_point_values"]
            # Alternative code for B1 sensitivity factor
            wmatK0 = zeros(ComplexF64, maxMeas, nPars + nNuisances)

            b1cp = 0.8:0.02:1.2    # set of b1 control points 
            for b1 in b1cp

                # simulate whole sequence for this B1 value
                # The line below deseveres some extra explanation.
                # In case of 'values', 
                # we will be calculating the discrepancy of the result by comparing m for B1 against m for B1=1.0,
                # and we assume the ρT1T2-Jacobian to be rather constant over this B1 region (but we take it halfway
                # between the two B1 values).
                b1midway = (B1metric == "multi_point_values") ? (1.0 + b1) / 2.0 : b1
                parameters = [BlochSimulators.T₁T₂B₁(T1test, T2test, b1midway)]
                m = simulate_magnetization(spgr, parameters)
                ∂m = simulate_derivatives_finite_difference(fit_parameters, m, spgr, parameters)
                wlocal[:, 1] = m
                wlocal[:, 2] = ∂m.T₁ .* T1test
                wlocal[:, 3] = ∂m.T₂ .* T2test
                wlocal[:, 4] = ∂m.B₁
                if (B1metric == "multi_point_values")
                    parameters_ideal = T₁T₂B₁(T1test, T2test, 1.0)
                    parameters_b1 = T₁T₂B₁(T1test, T2test, b1)
                    m_ideal = simulate_magnetization(spgr, parameters_ideal)
                    m_b1 = simulate_magnetization(spgr, parameters_b1)
                    wlocal[:, 4] = m_b1 .- m_ideal
                end

                # search for K=0 points an assemble W-mattrix for that set
                thisK0point = 0
                for i = 1:nTR
                    for sample in (trajectorySet[i])
                        if (0 <= sample.ky < 1) && (0 <= sample.kz < 1)
                            thisK0point += 1
                            if thisK0point <= maxMeas
                                wmatK0[thisK0point, :] = wlocal[i, :]
                            end
                        end
                    end
                end

                # calculate B1 sensitivity factor from W matrix
                W = wmatK0[:, 1:nPars]
                Wb1 = wmatK0[:, nPars+1]
                JhJ = adjoint(W) * W
                H0 = inv(JhJ + invReg[1:nPars, 1:nPars])
                b1fcpx = H0 * adjoint(W) * Wb1

                # accumulate effect of b1-sensitivities
                b1factors2 = b1factors2 .+ abs2.(b1fcpx)

            end # loop over B1 control points
            b1factors = sqrt.(b1factors2 ./ length(b1cp))
        else
            error("Unknown B1metric")
        end # condition on B1metric
    end # condition on b1 sensitivity analysis needed 

    return b1factors
end

function _analyze_jacobian(nky, nkz, nkyEff, useSym, nPars, nNuisances, nTR, wlocal, trajectorySet, options)

    maxMeas = options["maxMeas"]
    wmat = zeros(ComplexF64, nkyEff, nkz, maxMeas, nPars + nNuisances)

    lastMeas = zeros(Int64, nkyEff, nkz)
    for i = 1:nTR
        for sample in (trajectorySet[i])
            kyEff = sample.ky
            kzEff = sample.kz
            adjoin = (useSym && sample.ky < 0)
            if adjoin
                kzEff *= -1
                kyEff *= -1
            end

            floorky = Int(floor(kyEff))
            fracky = kyEff - floor(kyEff)
            wfracky = cos(fracky * pi / 2)
            floorkz = Int(floor(kzEff))
            frackz = kzEff - floor(kzEff)
            wfrackz = cos(frackz * pi / 2)

            for (thisky, wfracky) in [(floorky, cos(fracky * pi / 2)), (floorky + 1, sin(fracky * pi / 2))]
                for (thiskz, wfrackz) in [(floorkz, cos(frackz * pi / 2)), (floorkz + 1, sin(frackz * pi / 2))]
                    indexky = useSym ? thisky + 1 : thisky + nky ÷ 2 + 1
                    indexkz = thiskz + nkz ÷ 2 + 1
                    if ((indexky in 1:nkyEff) && (indexkz in 1:nkz))
                        thisMeas = lastMeas[indexky, indexkz] + 1
                        if thisMeas <= maxMeas
                            wmat[indexky, indexkz, thisMeas, :] = CondConj(adjoin, wlocal[i, :]) .* (wfracky * wfrackz)
                            lastMeas[indexky, indexkz] = thisMeas
                        end
                    end
                end
            end

            # enter the ky=kz=0 point twice, to select real part of noise only
            if (useSym && sample.ky == 0 && sample.kz == 0)
                thisky = 1
                thiskz = nkz ÷ 2 + 1
                thisMeas = lastMeas[thisky, thiskz] + 1
                if thisMeas <= maxMeas
                    wmat[thisky, thiskz, thisMeas, :] = conj.(wlocal[i, :])
                    lastMeas[thisky, thiskz] = thisMeas
                end
            end
        end
    end
    
    return wmat
end