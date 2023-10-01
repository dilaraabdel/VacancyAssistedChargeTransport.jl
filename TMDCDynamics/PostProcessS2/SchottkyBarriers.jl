

#=
Absolute change in Schottky barriers, absolute current and bias values w.r.t. time.
=#

module SchottkyBarriers

using ChargeTransport
using PyPlot
using DelimitedFiles
using ExtendableGrids

include("../TMDC.jl")

function main(;runSim = false, saveFig = false)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    params = "Params_TMDC_S2.jl"

    include("../../parameters/$params")

    if runSim
        sol, IV, biasVal = TMDC.main(;plotting = false, test = true, parameter_file = "../../parameters/$params", # "../../Params_TMDC_S2.jl"
                                        PrintBarriers = true)


        tvalues = sol.t; numbertsteps = length(tvalues)

        barrierChangeLeft =zeros(Float64, numbertsteps); barrierChangeRight = zeros(Float64, numbertsteps)

        for istep = 1:numbertsteps
            solution = sol[istep]
            ipsib    = 6           # this is the gradient of residual electric potential

            if solution[ipsib, 1] < 0
                barrierChangeLeft[istep]  = sqrt( - q/(4 * pi * ε0 * εi) * solution[ipsib, 1])
            else
                barrierChangeLeft[istep]  = 0.0
            end

            if solution[ipsib, end] < 0
                barrierChangeRight[istep] = sqrt( - q/(4 * pi * ε0 * εi) * solution[ipsib, end])
            else
                barrierChangeRight[istep] = 0.0
            end

        end

        writedlm("barrier-change-left-S2.dat",  [tvalues barrierChangeLeft])
        writedlm("barrier-change-right-S2.dat", [tvalues barrierChangeRight])
        writedlm("IV-S2-cycle-1-2-time.dat",    [tvalues IV])
        writedlm("bias-S2-cycle-1-2.dat",       [tvalues biasVal])

    else
        tvalues            = readdlm("data/barrier-change-left-S2.dat")[:, 1]
        barrierChangeLeft  = readdlm("data/barrier-change-left-S2.dat")[:, 2]
        barrierChangeRight = readdlm("data/barrier-change-right-S2.dat")[:, 2]
        IV                 = readdlm("data/IV-S2-cycle-1-2-time.dat")[:, 2]
        biasVal            = readdlm("data/bias-S2-cycle-1-2.dat")[:, 2]
    end




    PyPlot.semilogy(tvalues, abs.(Area .* IV), color = "darkblue", linewidth = 5)
    PyPlot.xlabel("time [s]", fontsize=17)
    PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)

    PyPlot.xlim(0.0, 16)
    PyPlot.ylim(1.0e-7, 2.0e-4)
    PyPlot.xticks([0, 5, 10, 15])
    PyPlot.yticks([1.0e-6, 1.0e-4])
    PyPlot.grid()
    PyPlot.tight_layout()
    if saveFig
        savefig("S2-time-current.pdf")
    end
    ##################################
    PyPlot.figure()
    PyPlot.plot(tvalues, barrierChangeLeft,  linewidth = 5, color="darkblue", label = "left")
    PyPlot.plot(tvalues, barrierChangeRight, linewidth = 5, color="dodgerblue", label = "right")

    PyPlot.xlabel("time [s]", fontsize=17)
    PyPlot.ylabel("\$ | \\Delta \\phi | \$[eV]", fontsize=17)
    PyPlot.xlim(0.0, 16)
    PyPlot.ylim(0.0, 0.14)
    PyPlot.xticks([0, 5, 10, 15])
    PyPlot.yticks([0.0, 0.1])
    PyPlot.grid()
    PyPlot.legend()
    PyPlot.tight_layout()

    if saveFig
        savefig("S2-time-barrier-change.pdf")
    end
    ##################################
    PyPlot.figure()
    PyPlot.plot(tvalues, biasVal,  linewidth = 5, color="gray", linestyle="--")

    PyPlot.xlabel("time [s]", fontsize=17)
    PyPlot.ylabel("bias [V]", fontsize=17)
    PyPlot.xlim(0.0, 16)
    PyPlot.ylim(-10, 10)
    PyPlot.xticks([0, 5, 10, 15])
    PyPlot.yticks([-10, 0, 10])
    PyPlot.grid()
    PyPlot.tight_layout()

    if saveFig
        savefig("S2-time-bias.pdf")
    end


end #  main


end # module
