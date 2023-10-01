
module Plotting_IV

#=
Comparison of I-V characteristics with measurements and for different cycles.

=#

using DelimitedFiles
using PyPlot

function main(saveFig = true)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    IV1     = readdlm("data/IV-S1-cycle-1.dat")
    IV2     = readdlm("data/IV-S1-cycle-2.dat")
    IV3     = readdlm("data/IV-S1-cycle-3.dat")
    IV4     = readdlm("data/IV-S1-cycle-4.dat")
    IV5     = readdlm("data/IV-S1-cycle-5.dat")
    IV6     = readdlm("data/IV-S1-cycle-6.dat")
    IV7     = readdlm("data/IV-S1-cycle-7.dat")
    IV8     = readdlm("data/IV-S1-cycle-8.dat")
    IV_meas = readdlm("../data/Li2018-FigS1i.dat")
    IV1NV   = readdlm("data/IV-S1-cycle-1-no-vacancies.dat")
    IV2NV   = readdlm("data/IV-S1-cycle-2-no-vacancies.dat")


    for xx in 1:length(IV_meas[:, 1])
        PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", markerfacecolor="none", markeredgewidth=2)
        if xx == length(IV_meas[:, 1])
            PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", label ="measured", markerfacecolor="none", markeredgewidth=2)
        end
    end
    PyPlot.semilogy(IV2[:, 1], IV2[:, 2], linewidth = 4, color = "darkblue",   label ="cycle 2")
    PyPlot.grid()
    PyPlot.legend(fontsize = 16)
    PyPlot.xlabel("bias [V]", fontsize=17)
    PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.xlim(-13, 13)
    PyPlot.ylim(1.0e-9, 1.0e-4)
    PyPlot.tight_layout()

    if saveFig
        savefig("S1-IV-measurement.pdf")
    end

    ##############################
    PyPlot.figure()
    PyPlot.semilogy(IV2[:, 1], IV2[:, 2], linewidth = 4, color = "darkblue",   label ="cycle 2-8")
    PyPlot.semilogy(IV3[:, 1], IV3[:, 2], linewidth = 4, color = "darkblue")
    PyPlot.semilogy(IV4[:, 1], IV4[:, 2], linewidth = 4, color = "darkblue")
    PyPlot.semilogy(IV5[:, 1], IV5[:, 2], linewidth = 4, color = "darkblue")
    PyPlot.semilogy(IV6[:, 1], IV6[:, 2], linewidth = 4, color = "darkblue")
    PyPlot.semilogy(IV7[:, 1], IV7[:, 2], linewidth = 4, color = "darkblue")
    PyPlot.semilogy(IV8[:, 1], IV8[:, 2], linewidth = 4, color = "darkblue")
    PyPlot.semilogy(IV1[:, 1], IV1[:, 2], linewidth = 4, color = "dodgerblue", label ="cycle 1")

    PyPlot.grid()
    PyPlot.legend(fontsize = 16)
    PyPlot.xlabel("bias [V]", fontsize=17)
    PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.xlim(-13, 13)
    PyPlot.ylim(1.0e-9, 1.0e-4)
    PyPlot.tight_layout()
    if saveFig
        savefig("S1-IV-cycle-comparison.pdf")
    end

    ##############################
    PyPlot.figure()

    for xx in 1:length(IV_meas[:, 1])
        PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", markerfacecolor="none", markeredgewidth=2)
        if xx == length(IV_meas[:, 1])
            PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", label ="measured", markerfacecolor="none", markeredgewidth=2)
        end
    end
    PyPlot.semilogy(IV1NV[:, 1], IV1NV[:, 2], linewidth = 4, color = "darkblue", label ="no vacancies (cycle 1-2)")
    PyPlot.semilogy(IV2NV[:, 1], IV2NV[:, 2], linewidth = 4, color = "darkblue")
    PyPlot.grid()
    PyPlot.legend(loc="lower left", fontsize = 16)
    PyPlot.xlabel("bias [V]", fontsize=17)
    PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()
    PyPlot.xlim(-13, 13)
    PyPlot.ylim(1.0e-9, 1.0e-4)

    if saveFig
        savefig("S1-IV-no-vacancies.pdf")
    end

    ##############################
end



end # module