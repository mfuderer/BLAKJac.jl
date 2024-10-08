# "conditional conjugation"
function CondConj(cond::Bool, f::Complex)
    cf = cond ? conj(f) : f
    return cf
end
function CondConj(cond::Bool, f::Vector{<:Complex})
    cf = cond ? conj.(f) : f
    return cf
end

function _plot_sequence(options, RFdeg, trajectorySet)

    plotOn = (length(options["plottypes"]) > 0)
    if plotOn
        options["plotfuncs"]["close"]()
    end
    if "first" ∈ options["plottypes"]
        options["plotfuncs"]["first"](RFdeg, trajectorySet, options)
    end
    if "trajectories" ∈ options["plottypes"]
        options["plotfuncs"]["trajectories"](trajectorySet)
    end
end

function _plot_weights(options, trajectorySet, wlocal)
    if "weights" ∈ options["plottypes"] # if "weights" in keys(options["plottypes"])
        kz = [(s)[1].kz for s in trajectorySet]
        ky = [(s)[1].ky for s in trajectorySet]
        #hue = 8*kz; 
        hue = ky
        options["plotfuncs"]["weighting"](wlocal[:, 1], wlocal[:, 2], wlocal[:, 3], hue, note)
    end
end

function _plot_original_jacobian(options, wlocal)
    if "original_jacobian" ∈ options["plottypes"]
        options["plotfuncs"]["originaljacobian"](wlocal[:, 1], wlocal[:, 2], wlocal[:, 3], note, options)
    end
end

function _plot_noise_spectrum(H, nes2, note, options)
    if "noiseSpectrum" ∈ options["plottypes"]
        options["plotfuncs"]["noisespectrum"](H, nes2, note, options)
    end
end

function _assemble_FISP3D(options)

    # TODO: This should be removed form BLAKJac: BLAKJac should work for any simulator and it is not the responsibility of BLAKJac to assemble the sequence)

    RFdeg = complex.(rand(options["nTR"])) 
    
    TR::Float64 = options["TR"]
    TE = TR / 2.01
    TI = (options["startstate"] == 1) ? 20.0 : 0.01

    cyclic::Bool = options["considerCyclic"]
    maxstate::Int64 = options["maxstate"]

    # additional parameter required for 3D simulations
    T_wait = 0.0 # 50.0# 0.0# 1.75 
    N_repeat = cyclic ? 5 : 1
    bINV = options["startstate"] < 0

    py_undersampling_factor = 1
    spgr = BlochSimulators.FISP3D(RFdeg, TR, TE, maxstate, TI, T_wait, N_repeat, bINV, false, py_undersampling_factor)
    return spgr
end

function _set_RF_train!(sequence, RFdeg)
    
    # Check whether RFdeg has the same element type (e.g. Complex) as what is expected by the sequence object
    if eltype(RFdeg) != eltype(sequence.RF_train)
        error("RFdeg and sequence.RF_train must have the same element type")
    end
    # Modify the sequence object's RF_train field in-place
    sequence.RF_train .= RFdeg

    return nothing
end

function _calculate_csf_penalty(options, sequence, T1T2set)

    if get(options, "lambda_CSF", 0.0) > 0.0
        
        # Set T₁ and T₂ values for CSF
        parameters = BlochSimulators.T₁T₂(4.0, 2.0)
        # Simulate the magnetization at echo times for CSF
        echos_csf = BlochSimulators.simulate_magnetization(cpu, sequence, [parameters])

        # Simulate for test values corresponding to different tissues as well (and sum?)
        echos_tissue = zero(echos_csf)
        for (index, (T1test, T2test)) in enumerate(T1T2set)
            parameters = BlochSimulators.T₁T₂(T1test, T2test)
            echos_tissue .+= BlochSimulators.simulate_magnetization(cpu, sequence, [parameters])
        end

        # Calculate the penalty
        CSF_penalty = norm(echos_csf) / (norm(echos_tissue) / length(T1T2set))

        # following vlock was intended for debugging, but the graph is nice enough to keep it
        if :PyPlot ∈ names(Main, imported=true)        
            i = options["optcount"]
            emergeCriterion = options["opt_emergeCriterion"] # (1*1000^3) ÷ (length(RFdegC)^2)
            emerge = (i % emergeCriterion == 2)
            if emerge
                Main.PyPlot.figure()
                Main.PyPlot.plot(abs.(echos_csf))
                for (index, (T1test, T2test)) in enumerate(T1T2set)
                    parameters = BlochSimulators.T₁T₂(T1test, T2test)
                    echos_tissue_one = BlochSimulators.simulate_magnetization(cpu, sequence, [parameters])
                    Main.PyPlot.plot(abs.(echos_tissue_one))
                end
                Main.PyPlot.pause(0.1)
            end
        end
        # end of debug part
    else
        CSF_penalty = 0.0
    end

    return CSF_penalty
end

function _plot_noise_levels(options, RFdeg, trajectorySet, noisesAll)

    RFrad = RFdeg * (π / 180)

    if haskey(options["plotfuncs"], "noisebars")
        options["plotfuncs"]["bars"](RFrad, trajectorySet, noisesAll)
    end
end

function _calculate_num_nuisance_parameters(options)
    
    useSurrogate = options["useSurrogate"]
    considerB1nuisance = false
    
    if (options["handleB1"] == "sensitivity") || (options["handleB1"] == "co-reconstruct")
        considerB1nuisance = true
    end
    nNuisances = (considerB1nuisance && !useSurrogate) ? 1 : 0
    
    return nNuisances
end