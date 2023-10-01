module PSC_comp_SteadyState

using PyPlot
using DelimitedFiles
using ChargeTransport
using ExtendableGrids

include("../parameters/Params_PSC_PCBM_MAPI_Pedot.jl")
include("PSCVolumeExclusion.jl")

# Eps set = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95]

function main(;saveFig = false, nonlinearDiff=true)

    tvaluesStead0p01, psil2Err0p01, iphinl2Err0p01, iphipl2Err0p01, iphial2Err0p01 = PSCVolumeExclusion.main(test = true, plotting = false,
                                                                                                    eps = 0.01, nonlinearDiffusion = nonlinearDiff,
                                                                                                    steadyStateCalc = true)

    tvaluesStead, psil2Err, iphinl2Err, iphipl2Err, iphial2Err = PSCVolumeExclusion.main(test = true, plotting = false,
                                                                eps = 0.7, nonlinearDiffusion = nonlinearDiff,
                                                                steadyStateCalc = true)


    function adjustWithTolerance!(xx)
        eps = 5.0e-18
        for i in eachindex(xx)
            if xx[i] <= eps
                xx[i] = eps
            end
        end
        return xx
    end

    plotStep = 11

    PyPlot.figure()
    subplot(211)

    ## 0p01
    PyPlot.semilogy(tvaluesStead0p01[1:plotStep:end], adjustWithTolerance!(psil2Err0p01[1:plotStep:end]),   color = "blue",  marker = "x", markersize=7,  linewidth =2, label = "\$ \\psi  - \\psi^\\infty \$")
    PyPlot.semilogy(tvaluesStead0p01[1:plotStep:end], adjustWithTolerance!(iphial2Err0p01[1:plotStep:end]), color = "gold", linestyle="dashed", marker = "o", linewidth =2, label = "\$ \\varphi_{\\mathrm{a}} - \\varphi_{\\mathrm{a}}^\\infty \$")

    PyPlot.ylabel("\$ || u - u^{\\infty} ||_{L^2}^2\$", fontsize=17)
    PyPlot.legend(fontsize=16)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.xlim(28, 130)
    PyPlot.ylim(1.0e-18, 1.0e-9)
    PyPlot.grid()
    PyPlot.tight_layout()
    #######################################################################
    #######################################################################
    subplot(212)

    ## 0p07
    PyPlot.semilogy(tvaluesStead[1:plotStep:end], adjustWithTolerance!(psil2Err[1:plotStep:end]),   color = "blue",  marker = "x", markersize=7, linewidth=2)
    PyPlot.semilogy(tvaluesStead[1:plotStep:end], adjustWithTolerance!(iphial2Err[1:plotStep:end]), color = "gold", linewidth=2, linestyle="dashed", marker = "o")

    PyPlot.xlabel("time [s]", fontsize=17)
    PyPlot.ylabel("\$ || u - u^{\\infty} ||_{L^2}^2\$", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.xlim(28, 130)
    PyPlot.ylim(1.0e-18, 1.0e-9)
    PyPlot.grid()
    PyPlot.tight_layout()

    if saveFig
        savefig("DD-parameter-nonlinearDiff-$nonlinearDiff-L2.pdf")
    end

end # main

end # module