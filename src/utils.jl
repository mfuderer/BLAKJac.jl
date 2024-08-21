function _plot_sequence(options, RFdeg, trajectorySet)

    plotOn = (length(options["plottypes"]) > 0)
    if plotOn
        options["plotfuncs"]["close"]()
    end
    if any(isequal("first"), options["plottypes"])
        options["plotfuncs"]["first"](RFdeg, trajectorySet, options)
    end
    if any(isequal("trajectories"), options["plottypes"])
        options["plotfuncs"]["trajectories"](trajectorySet)
    end
end

function _assemble_FISP3D(options, RFdeg)

    # TODO: This should be removed form BLAKJac: BLAKJac should work for any simulator and it is not the responsibility of BLAKJac to assemble the sequence)
    
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

    if any(i -> i == "noisebars", options["plottypes"])
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