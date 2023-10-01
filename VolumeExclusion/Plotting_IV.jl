

module Plotting_IV

using PyPlot
using DelimitedFiles

function main(;saveFig = false, withMeas = false)

    PyPlot.rc("font", family="serif", size=15)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    #########################
    if withMeas
        IV_meas     = readdlm("data/IV-measured-CaladoNature2016-forward-dark.dat")
        for ii in eachindex(IV_meas)
            if IV_meas[ii] < 0.0
                IV_meas[ii] = 0.0
            end
        end
    end

    IVjlMD0p01 = readdlm("data/0p01/drift/jl-IV-forward-drift-eps-0p01-no-generation.dat")
    IVjlMD0p5  = readdlm("data/0p5/drift/jl-IV-forward-drift-eps-0p5-no-generation.dat")
    IVjlMD0p9  = readdlm("data/0p9/drift/jl-IV-forward-drift-eps-0p9-no-generation.dat")

    IVjlND0p01 = readdlm("data/0p01/diffusion/jl-IV-forward-diffusion-eps-0p01-no-generation.dat")
    IVjlND0p5  = readdlm("data/0p5/diffusion/jl-IV-forward-diffusion-eps-0p5-no-generation.dat")
    IVjlND0p9  = readdlm("data/0p9/diffusion/jl-IV-forward-diffusion-eps-0p9-no-generation.dat")
    #########################

    IVimMD0p01 = readdlm("data/0p01/drift/IM-DD-parameter-J.dat")
    IVimMD0p5  = readdlm("data/0p5/drift/IM-DD-parameter-J.dat")
    IVimMD0p9  = readdlm("data/0p9/drift/IM-DD-parameter-J.dat")

    IVimND0p01 = readdlm("data/0p01/diffusion/IM-DD-parameter-J.dat")
    IVimND0p5  = readdlm("data/0p5/diffusion/IM-DD-parameter-J.dat")
    IVimND0p9  = readdlm("data/0p9/diffusion/IM-DD-parameter-J.dat")

    biasValues = IVjlMD0p01[:, 1]
    linewCT = 7; linewIM = 5

    if withMeas
        PyPlot.plot(IV_meas[:, 1], IV_meas[:, 2], linewidth = linewCT, color = "black", label = "measurement")
    end

    PyPlot.plot(biasValues, IVjlND0p01[:, 2]*(0.01^2).*1.0e3, linewidth = linewCT, marker = "o", color= "darkslategrey", label= "Nonlinear diffusion (\$\\epsilon = 0.01\$)")
    PyPlot.plot(biasValues, IVjlND0p5[:, 2]*(0.01^2).*1.0e3,  linewidth = linewCT, color= "mediumseagreen", label= "Nonlinear diffusion (\$\\epsilon = 0.5\$)")
    PyPlot.plot(biasValues, IVjlND0p9[:, 2]*(0.01^2).*1.0e3,  linewidth = linewCT, color= "mediumspringgreen", label= "Nonlinear diffusion (\$ \\epsilon = 0.9 \$)")

    PyPlot.plot(biasValues, -IVimND0p01[2:end], linestyle = ":", linewidth = linewIM, color = "black")
    PyPlot.plot(biasValues, -IVimND0p5[2:end],  linestyle = ":", linewidth = linewIM, color = "black")
    PyPlot.plot(biasValues, -IVimND0p9[2:end],  linestyle = ":", linewidth = linewIM, color = "black")

    PyPlot.grid()
    PyPlot.xlabel("bias [V]", fontsize=18)
    PyPlot.ylabel("current density [mA cm\$^{-2}\$ ]", fontsize=18)
    PyPlot.legend(loc = "upper left", fontsize=16)
    PyPlot.tick_params(axis="both", labelsize=18)

    if withMeas == false
        PyPlot.xlim(1.05, 1.205)
        PyPlot.ylim(-1.0, 35.0)
    else
        PyPlot.xlim(0.2, 1.23)
        PyPlot.ylim(-0.5, 7.0)
    end
    PyPlot.tight_layout()

    if saveFig
        if withMeas
            savefig("DD-IV-diffusion.pdf")
        else
            savefig("DD-IV-diffusion-zoom.pdf")
        end

    end

    ##########################################################################
    PyPlot.figure()

    if withMeas
        PyPlot.plot(IV_meas[:, 1], IV_meas[:, 2], linewidth = linewCT, color = "black", label = "measurement")
    end

    PyPlot.plot(biasValues, IVjlMD0p01[:, 2]*(0.01^2).*1.0e3, linewidth = linewCT, color= "maroon",    label= "Modified drift (\$\\epsilon = 0.01\$)")
    PyPlot.plot(biasValues, IVjlMD0p5[:, 2]*(0.01^2).*1.0e3,  linewidth = linewCT, color= "indianred", label= "Modified drift (\$\\epsilon = 0.5\$)")
    PyPlot.plot(biasValues, IVjlMD0p9[:, 2]*(0.01^2).*1.0e3,  linewidth = linewCT, color= "darksalmon", label= "Modified drift (\$ \\epsilon = 0.9 \$)")

    PyPlot.plot(biasValues, -IVimMD0p01[2:end], linestyle = ":", linewidth = linewIM, color = "black")
    PyPlot.plot(biasValues, -IVimMD0p5[2:end],  linestyle = ":", linewidth = linewIM, color = "black")
    PyPlot.plot(biasValues, -IVimMD0p9[2:end],  linestyle = ":", linewidth = linewIM, color = "black")

    PyPlot.grid()
    PyPlot.xlabel("bias [V]", fontsize=18)
    PyPlot.ylabel("current density [mA cm\$^{-2}\$ ]", fontsize=18)
    PyPlot.legend(loc = "upper left", fontsize=16)
    PyPlot.tick_params(axis="both", labelsize=18)

    if withMeas == false
        PyPlot.xlim(1.05, 1.205)
        PyPlot.ylim(-1.0, 35.0)
    else
        PyPlot.xlim(0.2, 1.23)
        PyPlot.ylim(-0.5, 7.0)
    end
    PyPlot.tight_layout()

    if saveFig
        if withMeas
            savefig("DD-IV-drift.pdf")
        else
            savefig("DD-IV-drift-zoom.pdf")
        end

    end

end

end
