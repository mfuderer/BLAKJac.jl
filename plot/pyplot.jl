function PlotClose()
    Main.PyPlot.close("all")
end

# plot of RF alongside with first elements of each trajectory
# intended for Cartesian measurements, neglecting kx ... so single-point
function PlotFirst(RFdeg::Vector{ComplexF64}, trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, options::Dict)
    # extract first element (presumably only element) of each trajectory
    ky = [s[1].ky for s in trajectorySet]
    kz = [s[1].kz for s in trajectorySet]
    nky = round(Int64, maximum(ky) - minimum(ky) + 1.0)
    nkz = round(Int64, maximum(kz) - minimum(kz) + 1.0)
    #figure(); plot(real.(RFdeg)); plot(ky);

    fig, ax1 = Main.PyPlot.plt.subplots()
    colorRF = [0, 0.2, 0.6]
    ax1.set_xlabel("TR index", fontsize=18)
    ax2 = ax1.twinx()
    ax2.set_ylabel("Flip angle (deg)", color=colorRF, fontsize=18)

    color = [1, 0.4, 0]
    pev = ky .- nky ÷ 2 .- 1.0
    ax1.set_ylabel("Phase encode values", color=color, fontsize=18)
    ax1.plot(pev, color=color)
    ax2.plot(real.(RFdeg), color=colorRF)
    #ax2.plot(abs.(RFdeg), color=colorRF)
    #ax2.plot(rad2deg.(angle.(RFdeg)), "--", color=colorRF)
    Main.PyPlot.title(options["rflabel"], fontsize=18)
    Main.PyPlot.plt.show()

    if nkz > 1
        color = [1, 0.2, 0]
        pev = kz .- nkz ÷ 2 .- 1.0
        ax1.plot(pev, color=color)
    end

    Main.PyPlot.fig.tight_layout()
    Main.PyPlot.plt.show()
end

# Show a set of (max e.g. 10) of the trajectories applied
function PlotTrajectories(trajectorySet::Vector{<:Vector{<:TrajectoryElement}})
    tl = length(trajectorySet)
    dl = min(tl, 10)
    Main.PyPlot.figure()
    for t in 1:dl
        trajectory = trajectorySet[t]
        ky = [s.ky for s in trajectory]
        kz = [s.kz for s in trajectory]
        Main.PyPlot.plot(ky, kz)
    end
end

# (from 2020)... plot sensitivity graphs
function PlotWeighing(w₀::Vector{ComplexF64}, w₁::Vector{ComplexF64}, w₂::Vector{ComplexF64}, hue::Vector{Float64}, note::String)
    σₘ = 1.0e-6
    invw₀ = [((abs(w) > σₘ) ? 1 / abs(w) : 1 / σₘ) for w in w₀]
    rs1 = abs.(w₁) .* invw₀
    rs2 = abs.(w₂) .* invw₀

    Main.PyPlot.figure()
    Main.PyPlot.plot(abs.(w₀ ./ w₀), label="relative sensitivity to rho")
    Main.PyPlot.plot(real.(rs1), label="relative sensitivity to T1")
    Main.PyPlot.plot(real.(rs2), label="relative sensitivity to T2")
    Main.PyPlot.title("$note")
    Main.PyPlot.legend()

    Main.PyPlot.figure()
    Main.PyPlot.plot(abs.(w₀), label="absolute sensitivity to rho")
    Main.PyPlot.plot(abs.(w₁), label="absolute sensitivity to T1")
    Main.PyPlot.plot(abs.(w₂), label="absolute sensitivity to T2")
    Main.PyPlot.title("$note")
    Main.PyPlot.legend()

    dotsize = 30.0 .* abs.(w₀)

    Main.PyPlot.figure()
    drs1 = min.(1.0, max.(0.0, real.(rs1)))
    drs2 = min.(1.0, max.(0.0, real.(rs2)))
    Main.PyPlot.scatter(drs1 .+ 0.0001 .* hue, drs2, dotsize, hue)
    Main.PyPlot.xlabel("relative sensitivity to T1 $note")
    Main.PyPlot.ylabel("relative sensitivity to T2")

    Main.PyPlot.figure()
    ds1 = min.(0.2, max.(-0.05, abs.(w₁)))
    ds2 = min.(0.2, max.(-0.05, abs.(w₂)))
    Main.PyPlot.scatter(ds1, ds2, dotsize, hue)
    Main.PyPlot.xlabel("absolute sensitivity to T1 $note")
    Main.PyPlot.ylabel("absolute sensitivity to T2")
end

function PlotOriginalJacobian(w₀::Vector{ComplexF64}, w₁::Vector{ComplexF64}, w₂::Vector{ComplexF64}, note::String, options)
    nky = 224
    jjj = zeros(length(w₀), 3)
    jjj[:, 1] = abs.(w₀)
    jjj[:, 2] = abs.(w₁)
    jjj[:, 3] = abs.(w₂)
    jjj = repeat(jjj, inner=(1, nky))
    Main.PyPlot.figure()
    Main.PyPlot.imshow(jjj)
end

function PlotBars(RFrad, trajectorySet::Vector{<:Vector{<:TrajectoryElement}}, noisesAll::Vector{Float64})
    ky = [(s)[1].ky for s in trajectorySet]
    rky = maximum(ky) - minimum(ky)
    rky = max(rky, 1.0)

    barnames = ["rho", "T1", "T2"]
    Main.PyPlot.figure()
    Main.PyPlot.xticks([1, 2, 3], barnames)
    Main.PyPlot.bar([1, 2, 3], noisesAll)
    Main.PyPlot.title("Normalized noise levels ")

    (fig, (ax1, ax2, ax3)) = (Main.PyPlot.subplots(1, 3, figsize=(8, 2)))

    ax3.set_xticks([1, 2, 3])
    ax3.set_xticklabels(barnames)
    ax3.bar([1, 2, 3], noisesAll)
    ax3.set_title("Normalized noise levels")

    ax1.plot(RFrad)
    ax1.plot(ky ./ rky)
    ax1.set(xlabel="nTR", ylabel="Flip angle (rad)")
    ax1.set_title("Simulated sequence")
end

function PlotNoiseSpectrum(H::Array{ComplexF64}, scaler, note::String, options)
    nPars = size(H, 3)
    nkz = size(H, 2)
    nkyEff = size(H, 1)
    kzprobe = 1 + nkz ÷ 2

    Hdiag = zeros(ComplexF64, nkyEff, nkz, nPars)
    for i in 1:nPars
        Hdiag[:, :, i] = H[:, :, i, i]
    end

    Main.PyPlot.figure()
    hrho = abs.(H[:, kzprobe, 1, 1])
    hT1 = abs.(H[:, kzprobe, 2, 2])
    hT2 = abs.(H[:, kzprobe, 3, 3])
    kyval = options["useSymmetry"] ? LinRange(1, nkyEff, nkyEff) : LinRange(-nkyEff / 2, nkyEff / 2, nkyEff)

    Main.PyPlot.plot(kyval, hrho, label="noise figure on rho")
    Main.PyPlot.plot(kyval, hT1, label="noise figure on T1")
    Main.PyPlot.plot(kyval, hT2, label="noise figure on T2")
    Main.PyPlot.legend(fontsize=16)
    Main.PyPlot.xlabel("spatial frequency", fontsize=16)
    Main.PyPlot.ylabel("noise variance (a.u.)", fontsize=16)
    Main.PyPlot.title(note)

    rootHmean = mean((sqrt.(abs.(Hdiag))), dims=3)
    scaledRootHmean = rootHmean[:, :, 1] .* sqrt(scaler)
    Main.PyPlot.figure()
    Main.PyPlot.imshow(scaledRootHmean)
    Main.PyPlot.title(note)
    @show mean(scaledRootHmean[2:nkyEff, :] .^ 2)
    @show mean(scaledRootHmean[1, :] .^ 2)
    @show mean(scaledRootHmean)
end

function PlotInfocon(I::Array{Float64}, note::String)
    Isum = sum(I, dims=3)
    if (size(Isum, 2) == 1)
        # no kz-extent; plot over ky
        Main.PyPlot.figure()
        Main.PyPlot.plot(Isum[:, 1])
        Main.PyPlot.title("Infocon map for $note")
        Main.PyPlot.figure()
        I_forShow = repeat(I[:, 1, :], inner=(1, 10))
        Main.PyPlot.imshow(I_forShow)
        Main.PyPlot.title("Infocon map for $note")
    else
        # show as image over (ky,kz)
        Main.PyPlot.figure()
        Main.PyPlot.imshow(Isum)
        Main.PyPlot.title("Infocon map for $note")
    end
end

function PlotIntermediate(RFdegC::Vector{ComplexF64}, options)
    if options["opt_complex"]
        Main.PyPlot.figure()
        Main.PyPlot.plot(abs.(RFdegC))
        Main.PyPlot.plot(57.3 * angle.(RFdegC))
    elseif options["opt_slow_phase"]
        fig, ax1 = Main.PyPlot.subplots()
        PlotAmplAndPhaseDeriv(RFdegC, ax1)
    else
        Main.PyPlot.figure()
        Main.PyPlot.plot(real.(RFdegC))
    end
    Main.PyPlot.pause(0.1)
end

function PlotAmplAndPhaseDeriv(RFdegC, ax1)
    # 'reconstruct' the second derivative of the phase
    d2 = similar(real.(RFdegC))
    d2[1:2] .= 0.0
    # division by 2 is to compensate for the square, which has been entered to remove phase=pi at 0-crossings
    d2[3:end] = [angle((RFdegC[i-2] * conj(RFdegC[i-1]^2) * RFdegC[i])^2) * 0.5 for i in 3:length(RFdegC)]

    color = "blue"
    ax1.set_xlabel("pulse number")
    ax1.set_ylabel("amplitude", color=color)
    ax1.plot(abs.(RFdegC), color=color)
    ax1.tick_params(axis="y", labelcolor=color)

    ax2 = ax1.twinx()  # instantiate a second Axes that shares the same x-axis

    color = "orange"
    ax2.set_ylabel("2nd derivative of phase", color=color)  # we already handled the x-label with ax1
    ax2.plot(rad2deg.(d2), color=color)
    ax2.tick_params(axis="y", labelcolor=color)
end
