function setPlotFuncs()

    if :PyPlot ∈ names(Main, imported=true)
        println("Use PyPlot-based plotting functions in BLAKJac")
        return Dict(
            "close" => PlotClose,
            "first" => PlotFirst,
            "trajectories" => PlotTrajectories,
            "weighing" => PlotWeighing,
            "originaljacobian" => PlotOriginalJacobian,
            "bars" => PlotBars,
            "noisespectrum" => PlotNoiseSpectrum,
            "infocon" => PlotInfocon,
            "intermediate" => PlotIntermediate
        )

    else
        println("No plotting in BLAKJac")
        return Dict(
            "close" => (x...) -> println("no plotting"),
            "first" => (x...) -> println("no plotting"),
            "trajectories" => (x...) -> println("no plotting"),
            "weighing" => (x...) -> println("no plotting"),
            "originaljacobian" => (x...) -> println("no plotting"),
            "bars" => (x...) -> println("no plotting"),
            "noisespectrum" => (x...) -> println("no plotting"),
            "infocon" => (x...) -> println("no plotting"),
            "intermediate" => (x...) -> println("no plotting")
        )
    end
end
