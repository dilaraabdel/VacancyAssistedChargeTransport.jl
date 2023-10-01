
module Plotting_IV

#=
Comparison of I-V characteristics with measurements and for different boundary contact cases.

=#

using DelimitedFiles
using PyPlot

function main(saveFig = true)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    IVSBL            = readdlm("data/IV-S2-cycle-2.dat")
    IVNoSBL          = readdlm("data/IV-S2-cycle-2-no-SBL.dat")
    IVNoSBLreduced   = readdlm("data/IV-S2-cycle-2-no-SBL-reduced-barrier.dat")
    ###################
    IVSBLNV          = readdlm("data/IV-S2-cycle-2-no-vacancies.dat")
    IVNoSBLNV        = readdlm("data/IV-S2-cycle-2-no-SBL-no-vacancies.dat")
    IVNoSBLreducedNV = readdlm("data/IV-S2-cycle-2-no-SBL-reduced-barrier-no-vacancies.dat")
    ###################
    IV_meas          = readdlm("../data/Li2018-FigS1j.dat")


    for xx in 1:length(IV_meas[:, 1])
        PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", markerfacecolor="none", markeredgewidth=2)
        if xx == length(IV_meas[:, 1])
            PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", label ="measured", markerfacecolor="none", markeredgewidth=2)
        end
    end
    PyPlot.semilogy(IVSBL[:, 1], IVSBL[:, 2], linewidth = 4, color = "darkblue",   label ="cycle 2")
    PyPlot.grid()
    PyPlot.legend(loc ="lower left", fontsize = 16)
    PyPlot.xlabel("bias [V]", fontsize=17)
    PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.xlim(-10, 10)
    PyPlot.ylim(1.0e-8, 5.0e-4)
    PyPlot.tight_layout()

    if saveFig
        savefig("S2-IV-measurement.pdf")
    end

    ##############################
    PyPlot.figure()
    for xx in 1:length(IV_meas[:, 1])
        PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", markerfacecolor="none", markeredgewidth=2)
        if xx == length(IV_meas[:, 1])
            PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", label ="measured", markerfacecolor="none", markeredgewidth=2)
        end
    end
    PyPlot.semilogy(IVNoSBL[:, 1], IVNoSBL[:, 2], linewidth = 4, color = "darkblue",   label ="without SBL")

    PyPlot.grid()
    PyPlot.legend(loc ="lower left", fontsize = 16)
    PyPlot.xlabel("bias [V]", fontsize=17)
    PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.xlim(-10, 10)
    PyPlot.ylim(1.0e-8, 5.0e-4)
    PyPlot.tight_layout()
    if saveFig
        savefig("S2-IV-without-SBL.pdf")
    end

    ##############################

    PyPlot.figure()

    for xx in 1:length(IV_meas[:, 1])
        PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", markerfacecolor="none", markeredgewidth=2)
        if xx == length(IV_meas[:, 1])
            PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "6", color = "gray", label ="measured", markerfacecolor="none", markeredgewidth=2)
        end
    end

    PyPlot.semilogy(IVSBLNV[:, 1],   IVSBLNV[:, 2],   linewidth = 4, color = "darkblue", label ="without vacancies")
    PyPlot.semilogy(IVNoSBLNV[:, 1], IVNoSBLNV[:, 2], linewidth = 4, color = "dodgerblue", linestyle = ":",  label ="without SBL & vacancies")

    PyPlot.grid()
    PyPlot.legend(loc ="lower left", fontsize = 16)
    PyPlot.xlabel("bias [V]", fontsize=17)
    PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.xlim(-10, 10)
    PyPlot.ylim(1.0e-8, 5.0e-4)
    PyPlot.tight_layout()

    if saveFig
        savefig("S2-IV-no-vacancies.pdf")
    end

    ##############################

    PyPlot.figure()

    PyPlot.semilogy(IVSBL[:, 1], IVSBL[:, 2], linewidth = 4, color = "darkblue",   label ="with SBL")
    PyPlot.semilogy(IVNoSBLreduced[:, 1], IVNoSBLreduced[:, 2], linewidth = 4, color = "dodgerblue", linestyle = "--",  label ="no SBL, reduced barrier")

    PyPlot.grid()
    PyPlot.legend(loc ="lower left", fontsize = 16)
    PyPlot.xlabel("bias [V]", fontsize=17)
    PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.xlim(-10, 10)
    PyPlot.ylim(1.0e-8, 5.0e-4)
    PyPlot.tight_layout()

    if saveFig
        savefig("S2-IV-reduced-barrier.pdf")
    end

end



end # module