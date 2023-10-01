
module Plotting_BandEdges

#=
Comparison of band edges for the first and the second cycle.

=#

using DelimitedFiles
using PyPlot
using ChargeTransport

include("../../parameters/Params_TMDC_S1.jl")

function main(saveFig = true)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    sol1 = readdlm("data/sol-S1-cycle-1.dat")
    sol2 = readdlm("data/sol-S1-cycle-2.dat")
    sol3 = readdlm("data/sol-S1-cycle-3.dat")

    coord = sol1[:, 1]
    iphin = 2
    iphip = 3
    iphia = 4
    ipsi  = 5

    PyPlot.plot(coord./μm, En/q .- sol1[:, ipsi], color = "green", linewidth = 4, label = " \$ E_{\\mathrm{n}} - q\\psi \$")
    PyPlot.plot(coord./μm, Ep/q .- sol1[:, ipsi], color = "red",   linewidth = 4, label = " \$ E_{\\mathrm{p}} - q\\psi \$")
    PyPlot.plot(coord./μm, -Ea/q .+ sol1[:, ipsi], color = "gold", linewidth = 4, label = " \$ -E_{\\mathrm{a}} + q\\psi \$")
    PyPlot.plot(coord./μm, -sol1[:, iphin], color = "black",       linewidth = 4, label = " \$ -q \\varphi_{\\mathrm{n}, \\mathrm{p}}\$", linestyle ="--")
    PyPlot.plot(coord./μm,  sol1[:, iphia], color = "gold",       linewidth = 5, label = " \$ q \\varphi_{\\mathrm{a}}\$", linestyle =":")

    PyPlot.grid()
    PyPlot.legend(fontsize=14)
    PyPlot.xticks( [0.0, 0.5, 1.0] )
    PyPlot.yticks( [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5] )
    PyPlot.ylim(-1.5, 0.5)
    PyPlot.title("\$ t = 0 \$ s, \$ U = 0 \$ V")
    PyPlot.xlabel("length [\$ \\mu\$m]", fontsize=17)
    PyPlot.ylabel("energy [eV]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig("S1-banddiagram-cycle-1.pdf")
    end

    ###################################
    PyPlot.figure()

    PyPlot.plot(coord./μm, En/q .- sol2[:, ipsi], color = "green", linewidth = 4, label = " \$ E_{\\mathrm{n}} - q\\psi \$")
    PyPlot.plot(coord./μm, Ep/q .- sol2[:, ipsi], color = "red",   linewidth = 4, label = " \$ E_{\\mathrm{p}} - q\\psi \$")
    PyPlot.plot(coord./μm, -Ea/q .+ sol2[:, ipsi], color = "gold", linewidth = 4, label = " \$ -E_{\\mathrm{a}} + q\\psi \$")
    PyPlot.plot(coord./μm, -sol2[:, iphin], color = "black",       linewidth = 4, label = " \$ -q \\varphi_{\\mathrm{n}, \\mathrm{p}}\$", linestyle ="--")
    PyPlot.plot(coord./μm,  sol2[:, iphia], color = "gold",       linewidth = 5, label = " \$ q \\varphi_{\\mathrm{a}}\$", linestyle =":")

    PyPlot.grid()
    PyPlot.legend(fontsize=14)
    PyPlot.yticks( [-1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5] )
    PyPlot.ylim(-1.5, 0.5)
    PyPlot.xticks( [0.0, 0.5, 1.0] )
    PyPlot.title("\$ t = 10.4 \$ s, \$ U = 0 \$ V")
    PyPlot.xlabel("length [\$ \\mu\$m]", fontsize=17)
    PyPlot.ylabel("energy [eV]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()

    if saveFig
        savefig("S1-banddiagram-cycle-2.pdf")
    end

    ###################################
    # PyPlot.figure()

    # PyPlot.plot(coord./μm, En/q .- sol3[:, ipsi], color = "green", linewidth = 4, label = " \$ E_{\\mathrm{n}} - q\\psi \$")
    # PyPlot.plot(coord./μm, Ep/q .- sol3[:, ipsi], color = "red",   linewidth = 4, label = " \$ E_{\\mathrm{p}} - q\\psi \$")
    # PyPlot.plot(coord./μm, -Ea/q .+ sol3[:, ipsi], color = "gold", linewidth = 4, label = " \$ -E_{\\mathrm{a}} + q\\psi \$")
    # PyPlot.plot(coord./μm, -sol3[:, iphin], color = "black",       linewidth = 4, label = " \$ -q \\varphi_{\\mathrm{n}, \\mathrm{p}}\$", linestyle ="--")
    # PyPlot.plot(coord./μm,  sol3[:, iphia], color = "gold",       linewidth = 4, label = " \$ q \\varphi_{\\mathrm{a}}\$", linestyle =":")

    # PyPlot.grid()
    # PyPlot.legend(fontsize=14)
    # PyPlot.title("\$ t = 20.8 \$ s, \$ U = 0 \$ V")
    # PyPlot.xlabel("length [\$ \\mu\$m]", fontsize=17)
    # PyPlot.ylabel("energy [eV]", fontsize=17)
    # PyPlot.tick_params(which ="both", labelsize=18)
    # PyPlot.tight_layout()
end



end # module