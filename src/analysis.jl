
"""
BLAKJac_analysis!(resource::CPU1, RFdeg::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict, saved_H::Dict=Dict())

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

function BLAKJac_analysis!(resource::CPU1, RFdeg::Vector{ComplexF64}, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict, saved_H::Dict=Dict())

    RFrad = RFdeg * (π / 180)
    TR::Float64 = options["TR"]
    TE = TR / 2.01
    TI = (options["startstate"] == 1) ? 20.0 : 0.01
    T1ref::Float64 = options["T1ref"]
    T2ref::Float64 = options["T2ref"]
    maxstate::Int64 = options["maxstate"]
    useSurrogate::Bool = options["useSurrogate"]
    T1T2set::Vector{Tuple{Float64,Float64}} = options["T1T2set"]
    cyclic::Bool = options["considerCyclic"]

    nPars = 3
    considerB1nuisance = false
    if (options["handleB1"] == "sensitivity") || (options["handleB1"] == "co-reconstruct")
        considerB1nuisance = true
    end
    nNuisances = (considerB1nuisance && !useSurrogate) ? 1 : 0

    # plot simulated sequence if requested
    plotOn = (length(options["plottypes"]) > 0)
    if plotOn
        options["plotfuncs"]["close"]()
    end
    if any(i -> i == "first", options["plottypes"])
        options["plotfuncs"]["first"](RFdeg, trajectorySet, options)
    end
    if any(i -> i == "trajectories", options["plottypes"])
        options["plotfuncs"]["trajectories"](trajectorySet)
    end

    # Assemble SPGR sequence simulator
    # RFdegExtended = repeat(RFdeg,outer=(cyclic ? 2 : 1))

    # additional parameter required for 3D simulations
    T_wait = 0.0 # 50.0# 0.0# 1.75 
    N_repeat = cyclic ? 5 : 1
    bINV = options["startstate"] < 0
    py_undersampling_factor = 1
    spgr = BlochSimulators.FISP3D(RFdeg, TR, TE, maxstate, TI, T_wait, N_repeat, bINV, false, py_undersampling_factor)

    # Initialize accumulators to calculate averaged result over the T1T2set
    b1factors2All = zeros(nPars)
    noisesAll = zeros(nPars)
    b1factorsLast = zeros(nPars)
    ItotAll = 0.0

    # Loop over all (e.g. 7) probe-values of (T1,T2)
    for (index, (T1test, T2test)) in enumerate(T1T2set)
        B1test = 1.0

        H, noises, Itot, b1factors = BLAKJacOnSingleT1T2(resource, T1test, T2test, B1test, nNuisances, spgr, trajectorySet, options)

        # accumulate
        b1factors2All .+= (b1factors) .^ 2
        b1factorsLast = b1factors
        noisesAll .+= noises
        ItotAll += Itot

        if (@isdefined saved_H)
            saved_H["$index for $(options["rfFile"])"] = H
        end
    end #for loop over probe-set of (T1,T2)

    # ------------------------------------- estimate of value for CSF_penalty
    if get(options, "lambda_CSF", 0.0) > 0.0
        T₁T₂_csf = T₁T₂(4.0, 2.0)
        echos_csf = BlochSimulators.simulate_magnetization(spgr, T₁T₂_csf)

        echos_tissue = zeros(ComplexF64, size(echos_csf))
        for (index, (T1test, T2test)) in enumerate(T1T2set)
            T₁T₂_tissue = T₁T₂(T1test, T2test)
            echos_tissue .+= BlochSimulators.simulate_magnetization(spgr, T₁T₂_tissue)
        end
        CSF_penalty = norm(echos_csf) / (norm(echos_tissue) / length(T1T2set))

        # following vlock was intended for debugging, but the graph is nice enough to keep it        
        i = options["optcount"]
        emergeCriterion = options["opt_emergeCriterion"] # (1*1000^3) ÷ (length(RFdegC)^2)
        emerge = (i % emergeCriterion == 2)
        if emerge
            Main.PyPlot.figure()
            Main.PyPlot.plot(abs.(echos_csf))
            for (index, (T1test, T2test)) in enumerate(T1T2set)
                T₁T₂_tissue_one = T₁T₂(T1test, T2test)
                echos_tissue_one = BlochSimulators.simulate_magnetization(sequence, T₁T₂_tissue_one)
                Main.PyPlot.plot(abs.(echos_tissue_one))
            end
            Main.PyPlot.pause(0.1)
        end
        # end of debug part


    else
        CSF_penalty = 0.0
    end

    b1factorsRms = length(T1T2set) == 1 ? b1factorsLast : sqrt.(b1factors2All ./ size(T1T2set))
    noisesAll ./= size(T1T2set)
    ItotAll /= only(size(T1T2set))

    if any(i -> i == "noisebars", options["plottypes"])
        options["plotfuncs"]["bars"](RFrad, trajectorySet, noisesAll)
    end

    return noisesAll, ItotAll, b1factorsRms, CSF_penalty
end

function BLAKJacOnSingleT1T2(resource::CPU1, T1test, T2test, B1test, nNuisances, spgr::BlochSimulators.FISP3D, trajectorySet::Vector{Vector{TrajectoryElement}}, options::Dict)
    TR = options["TR"]
    nTR = length(trajectorySet)
    nky = options["nky"]
    nkz = options["nkz"]
    useSym = options["useSymmetry"]
    nkyEff = useSym ? nky ÷ 2 : nky
    nkyEff = max(1, nkyEff)
    nPars = 3

    T1ref::Float64 = options["T1ref"]
    T2ref::Float64 = options["T2ref"]
    useSurrogate = options["useSurrogate"]
    cyclic = options["considerCyclic"]
    note = @sprintf("for (T1,T2)=(%6.1f,%6.1f)", 1000 * T1test, 1000 * T2test)
    takeB1asVariable = (options["handleB1"] == "co-reconstruct")

    parameters = (nNuisances > 0) ? T₁T₂B₁(T1test, T2test, B1test) : T₁T₂(T1test, T2test)
    parameters = StructVector([parameters])
    fit_parameters = (nNuisances > 0) ? (:T₁, :T₂, :B₁) : (:T₁, :T₂)

    # Now an inelegant if-then-else approach follows,
    # prompted by not fully understanding the struct that is returned by simulate_derivatives() ...
    wlocal = zeros(ComplexF64, nTR, nPars + nNuisances)

    if (useSurrogate)
        error("Surrogate model not supported at the moment")
        # surrogate = BlochSimulators.PolynomialSurrogate(spgr, options)
        # derivs = simulate_derivatives(resource, surrogate, parameters, fit_parameters)
        # wlocal[:, 1] = derivs.m
        # wlocal[:, 2] = derivs.∂T₁
        # wlocal[:, 3] = derivs.∂T₂
    else

        m = simulate_magnetization(spgr, parameters)
        ∂m = simulate_derivatives_finite_difference(fit_parameters, m, spgr, parameters)

        wlocal[:, 1] = m
        wlocal[:, 2] = ∂m.T₁ .* T1test
        wlocal[:, 3] = ∂m.T₂ .* T2test
        if nNuisances > 0
            wlocal[:, 4] = ∂m.B₁
        end
    end

    if any(i -> i == "weights", options["plottypes"]) # if "weights" in keys(options["plottypes"])
        kz = [(s)[1].kz for s in trajectorySet]
        ky = [(s)[1].ky for s in trajectorySet]
        #hue = 8*kz; 
        hue = ky
        options["plotfuncs"]["weighting"](wlocal[:, 1], wlocal[:, 2], wlocal[:, 3], hue, note)
    end

    if any(i -> i == "original_jacobian", options["plottypes"])
        options["plotfuncs"]["originaljacobian"](wlocal[:, 1], wlocal[:, 2], wlocal[:, 3], note, options)
    end
    # plotFuncs["jacobian"]()

    # analyze Jacobian
    invRegB1 = takeB1asVariable ? 0.0 : sqrt(prevfloat(Inf)) # corretion 2024-06-04. Logic: if B1 is to be co-reconstructed, 
    # then it is not regularized; if it is to be sensitivity-analyzed, 
    # it should not influence the estimate of T1T2rho noise, so it is
    # very heavily ('infinitely') regularized
    invReg = Diagonal([options["invregval"]; [invRegB1]])
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

    nParsX = nPars
    nNuisancesX = nNuisances
    if (nPars == 3) && (nNuisances == 1) && takeB1asVariable
        nParsX = 4
        nNuisancesX = 0
    end

    # analyze matrices for all ky
    sumH = zeros(ComplexF64, nParsX, nParsX)
    H = zeros(ComplexF64, nkyEff, nkz, nParsX, nParsX)
    Hdiag = zeros(ComplexF64, nkyEff, nkz, nParsX)
    I = zeros(Float64, nkyEff, nkz, nParsX)
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

    # B1 sensitivity factor
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

    sumH = sum(H_emphasized, dims=1)
    sumH = sum(sumH, dims=2)

    # A scaling factor is introduced that should normalize the noise level to 1 in case of Rho-only reconstruction
    # "nes2"= Normalized Expected Signal squared
    # Still to be added to notebook per 2021-02-01
    nes2 = 2.0 * 4.0 * (nTR * TR + T1ref) * T2ref / (nky * nkz * (T1ref + T2ref)^2) # questionable whether ok if useSym
    # The extra factor of 2.0, is yet to be explained

    for i in 1:nPars
        Hdiag[:, :, i] = H[:, :, i, i]
    end

    # output diagonal elements
    if any(i -> i == "noiseSpectrum", options["plottypes"])
        options["plotfuncs"]["noisespectrum"](H, nes2, note, options)
    end

    # output RMS of Hessian
    sigmaRho = sqrt(abs(sumH[1, 1, 1, 1]) / nkyEff / nkz)
    sigmaT1 = sqrt(abs(sumH[1, 1, 2, 2]) / nkyEff / nkz)
    sigmaT2 = sqrt(abs(sumH[1, 1, 3, 3]) / nkyEff / nkz)
    noises = [sigmaRho, sigmaT1, sigmaT2]
    noises .*= sqrt(nes2)

    # info content
    σref2 = options["sigma_ref"]^2

    kym = useSym ? 1 : nky ÷ 2 + 1
    kzm = nkz ÷ 2 + 1
    pevaly = [((i - kym)^2 + 1^2) / ((nky ÷ 2)^2 + 1^2) for i in 1:nkyEff]
    pevalz = [((i - kzm)^2 + 1^2) / ((nkz ÷ 2)^2 + 1^2) for i in 1:nkz]

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
    if any(i -> i == "infocon", options["plottypes"])
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

    return H, noises, Itot, b1factors
end
