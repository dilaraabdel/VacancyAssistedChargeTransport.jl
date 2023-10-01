
module Plotting_BandEdges

#=
Comparison of band edges for the first cycle with and without SBL.

=#

using DelimitedFiles
using PyPlot
using ChargeTransport

include("../../parameters/Params_TMDC_S2.jl")

function main(saveFig = true)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    solSBL   = readdlm("data/sol-S2-cycle-1.dat")
    solNoSBL = readdlm("data/sol-S2-cycle-1-no-SBL.dat")

    coord = solSBL[:, 1]
    ipsi  = 5

    subplot(211)

    PyPlot.plot(coord./μm, En/q .- solNoSBL[:, ipsi], color = "dodgerblue", linewidth = 5, label = " without SBL")
    PyPlot.plot(coord./μm, En/q .- solSBL[:, ipsi], color = "darkblue", linewidth = 5, label = " with SBL")

    PyPlot.grid()
    PyPlot.legend(fontsize=14)
    PyPlot.xlim(-0.01, 0.025)
    PyPlot.xticks( [0.0, 0.025] )
    PyPlot.ylim(-0.0, 0.2)
    PyPlot.title("\$ t = 0 \$ s, \$ U = 0 \$ V")
    PyPlot.xlabel("length [\$ \\mu\$m]", fontsize=17)
    PyPlot.ylabel("\$ E_{\\mathrm{n}} - q\\psi \$ [eV]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    subplot(212)

    PyPlot.plot(coord./μm, En/q .- solNoSBL[:, ipsi], color = "dodgerblue", linewidth = 5, label = " without SBL")
    PyPlot.plot(coord./μm, En/q .- solSBL[:, ipsi], color = "darkblue", linewidth = 5, label = " with SBL")

    PyPlot.grid()
    PyPlot.legend(fontsize=14)
    PyPlot.xlim(0.975, 1.01)
    PyPlot.xticks( [0.975, 1.0] )
    PyPlot.ylim(0.0, 0.2)
    PyPlot.xlabel("length [\$ \\mu\$m]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()
    if saveFig
        savefig("S2-banddiagram-cycle-1.pdf")
    end

end



end # module