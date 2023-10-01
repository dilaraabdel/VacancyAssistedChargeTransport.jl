

#=
Comparison of I-V curves for different mobilities
=#

module DifferentMobilityIV

using ChargeTransport
using PyPlot
using DelimitedFiles
using ExtendableGrids

include("../TMDC.jl")

function main(;runSim = false)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    mobilities = [1.0e-15, 1.0e-14, 3.0e-14, 5.0e-14, 9.0e-14, 2.0e-13] * (m^2) / (V*s)


    for μa in mobilities

        helpμa   = collect(string(μa)); helpμa[2] = 'p'
        textμa   = join(helpμa)

        println("         μa = $μa")

        if runSim
            biasVal, IV = TMDC.main(;plotting = false, test = true, parameter_file = "../../parameters/Params_TMDC_S1.jl", # "../../parameters/Params_TMDC_S2.jl"
                                    mobilitySweep           = true,  μaSweep       = μa)

            writedlm("IV-S1-cycle-2-mobility-$textμa.dat", [biasValues IV])
        else
            data    = readdlm("data/IV-S1-cycle-2-mobility-$textμa.dat")
            biasVal = data[:, 1]
            IV      = data[:, 2]
        end

        PyPlot.figure()
        PyPlot.semilogy(biasVal, IV, linewidth = 7, color = "dodgerblue")
        PyPlot.grid()
        PyPlot.title("\$ \\mu_{\\mathrm{a}} = \$ $μa \$  \\mathrm{m}^2/\$(Vs)")
        PyPlot.xlabel("bias [V]", fontsize=17)
        PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.xlim(-13, 13)
        PyPlot.ylim(1.0e-8, 1.0e-4)
        PyPlot.xticks([-10, -5, 0, 5, 10])
        PyPlot.yticks([1.0e-8, 1.0e-6, 1.0e-4])
        PyPlot.tight_layout()
        PyPlot.gcf()
        savefig("IV-S1-cycle-2-mobility-$textμa.pdf")
    end



end #  main


end # module
