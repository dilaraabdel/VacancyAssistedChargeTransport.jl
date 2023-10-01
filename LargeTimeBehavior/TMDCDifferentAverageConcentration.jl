

#=

=#

module TMDCDifferentAverageConcentration

using ChargeTransport
using PyPlot
using DelimitedFiles
using ExtendableGrids

include("TMDC.jl")
include("../parameters/Params_TMDC_S1.jl")

function main(;test = false, runSim = false, saveFig = false)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    ## vector containing different values for Ea resulting in different average vacancy concentrations
    EaLoop = collect(-5.02:0.005:-3.37).* eV

    if runSim

        trelEntropVec = zeros(0); tpsiErrVec  = zeros(0); tphinErrVec = zeros(0);
        tphipErrVec   = zeros(0); tphiaErrVec = zeros(0); naAvgVec    = zeros(0)

        for Ea in EaLoop

            if test == false
                println("         Ea = $(Ea/q) eV")
            end

            try
                tval, naAvg = TMDC.main(;plotting   = false, test    = true,
                                         naAvgSweep = true,  EaSweep = Ea)

                push!(trelEntropVec, tval[1]); push!(tpsiErrVec, tval[2]); push!(tphinErrVec, tval[3]);
                push!(tphipErrVec, tval[4]); push!(tphiaErrVec, tval[5]); push!(naAvgVec, naAvg);
            catch y
                @warn("Exception: ", y) # only warn and continue in case of an error
            end

        end

        writedlm("data/naAvg.dat",       naAvgVec);    writedlm("data/trelEntropVec.dat", trelEntropVec)
        writedlm("data/tpsiErrVec.dat",  tpsiErrVec);  writedlm("data/tphinErrVec.dat",   tphinErrVec)
        writedlm("data/tphipErrVec.dat", tphipErrVec); writedlm("data/tphiaErrVec.dat",   tphiaErrVec)

    else
        naAvgVec    = readdlm("data/naAvg.dat");       trelEntropVec = readdlm("data/trelEntropVec.dat")
        tpsiErrVec  = readdlm("data/tpsiErrVec.dat");  tphinErrVec   = readdlm("data/tphinErrVec.dat")
        tphipErrVec = readdlm("data/tphipErrVec.dat"); tphiaErrVec   = readdlm("data/tphiaErrVec.dat")
    end

    plotStep = 1:4:length(naAvgVec)

    PyPlot.semilogx(naAvgVec[plotStep], trelEntropVec[plotStep],  color = "darkgreen", marker = "x", label = "relative entropy w.r.t. steady state")
    PyPlot.xlabel("average vacancy concentration [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
    PyPlot.ylabel("convergence time [s]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.grid()
    PyPlot.legend(loc = "upper left")
    PyPlot.ylim(0.0, 220)
    PyPlot.tight_layout()
    if saveFig
        savefig("rel-entropy-TMDC-different-naAvg.pdf")
    end

    #####################################################
    PyPlot.figure()

    PyPlot.semilogx(naAvgVec[plotStep], tpsiErrVec[plotStep], marker = "o", markersize=11, color = "blue", linewidth =4, label = "\$ || \\psi  - \\psi^\\infty ||_{L^2}^2 \$")
    PyPlot.semilogx(naAvgVec[plotStep], tphinErrVec[plotStep], marker = "o", markersize=10,color = "green",  linewidth =4,  label = "\$||  \\varphi_{\\mathrm{n}} - \\varphi_{\\mathrm{n}}^\\infty  ||_{L^2}^2\$")
    PyPlot.semilogx(naAvgVec[plotStep], tphipErrVec[plotStep], marker = "o", markersize=7,color = "red", linewidth =4,  linestyle=":",  label = "\$ || \\varphi_{\\mathrm{p}} - \\varphi_{\\mathrm{p}}^\\infty ||_{L^2}^2 \$")
    PyPlot.semilogx(naAvgVec[plotStep], tphiaErrVec[plotStep], marker = "o", markersize=4,color = "gold", label = "\$ || \\varphi_{\\mathrm{a}} - \\varphi_{\\mathrm{a}}^\\infty ||_{L^2}^2 \$")

    PyPlot.xlabel("average vacancy concentration [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
    PyPlot.ylabel("convergence time [s]", fontsize=17)
    PyPlot.legend(fontsize=16)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.grid()
    PyPlot.ylim(0.0, 220)
    PyPlot.tight_layout()
    if saveFig
        savefig("l2-error-TMDC-different-naAvg.pdf")
    end

end #  main


end # module
