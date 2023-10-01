

#=
Calculate steady state solution for a given threshold voltage
=#

module SteadyState

using ChargeTransport
using PyPlot

include("../TMDC.jl")

function main(;saveFig = true)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Calculate quasi Static solution")
    end
    ################################################################################

    biasValSteadyStat = 9.6 * V # 1.6 * V # 5.1 * V, 11.2 * V, 9.6 * V, 4.4 * V, -2.6 * V
    coord, sol, data = TMDC.main(;plotting = false, test = true, parameter_file = "../../parameters/Params_TMDC_S1.jl", # "../../parameters/Params_TMDC_S2.jl"
                                  SteadyStatCalculation = true, biasValSteadyStat = biasValSteadyStat)


    ################################################################################
    if test == false
        println("Some plotting")
    end
    ################################################################################

    iphia       = 3
    regionflake = 1
    inode       = 1:length(coord)

    inode       = 1:length(coord)
    na          = get_density(sol[1],   data, iphia, regionflake, inode = inode)
    naStead     = get_density(sol[end], data, iphia, regionflake, inode = inode)

    PyPlot.semilogy(coord./μm, naStead, color = "silver", linewidth = 7,  label = "steady state solution")
    PyPlot.semilogy(coord./μm, na,      color = "dodgerblue", linewidth = 7,  label = " \$U\$ = $biasValSteadyStat V")

    PyPlot.grid()
    PyPlot.xlabel("length [\$ \\mu \$m]", fontsize=17)
    PyPlot.ylabel("vacancy density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize=16)
    PyPlot.xlim(0.0, 1.0)
    PyPlot.ylim(1.0e15, 1.0e25)
    PyPlot.xticks([0.0, 0.5, 1.0])
    PyPlot.yticks([1.0e15, 1.0e20, 1.0e25])
    PyPlot.tight_layout()
    if saveFig
        savefig("steady-state-$biasValSteadyStat-V.pdf")
    end

    testval = sum(filter(!isnan, sol[end]))/length(sol[end]) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
      testval = -1.106440770523787
    main(test = true) ≈ testval
end


end # module
