
module PSC_comp_one_eps

using PyPlot
using DelimitedFiles
using ChargeTransport
using ExtendableGrids

include("../parameters/Params_PSC_PCBM_MAPI_Pedot.jl")
include("PSCVolumeExclusion.jl")

# Eps set = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95]

function main(;saveFig = false, eps = 0.01)

    PyPlot.rc("font", family="serif", size=15)
    PyPlot.rc("mathtext", fontset="dejavuserif")

    helpEps   = collect(string(eps)); helpEps[2] = 'p'
    textEps   = join(helpEps)

    ## read in the grid from Ionmonger
    shift     = readdlm("data/IM-DD-parameter-grid-ETL.dat")[1]
    grid_ETL  = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-ETL.dat") .- shift)
    length_n  = length(grid_ETL)
    grid_intr = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-intr.dat").- shift)
    length_ni = length_n + length(grid_intr)
    grid_HTL  = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-HTL.dat") .- shift)

    coord     = glue(grid_ETL[:, 1], grid_intr[:, 1])
    coord     = glue(coord[:, 1], grid_HTL[:, 1])

    ##############################################################################################
    ############################       nonlinear diffusion data       ############################
    ##############################################################################################

    folder           = "diffusion"

    ## built-in voltage
    UT               = kB * T / q
    Vbi              = (En[regionDonor] - Ep[regionAcceptor] - kB * T * log(Nn[regionDonor] * Np[regionAcceptor]/ (Cn * Cp)))/q
    Nintr_d          = sqrt(Nn[regionDonor] * Np[regionDonor]) * exp(-(En[regionDonor] - Ep[regionDonor] ) / (2 * kB * T))
    psi0_d           = (En[regionDonor]+ Ep[regionDonor])/(2 * q) - 0.5 * UT * log(Nn[regionDonor]/Np[regionDonor]) + UT * asinh(Cn/(2*Nintr_d))
    shiftPsi         = - Vbi/2 + psi0_d

    psiIM_ETL0p0ND   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-0p0.dat")   .+ shiftPsi
    psiIM_INTR0p0ND  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-0p0.dat")  .+ shiftPsi
    psiIM_HTL0p0ND   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-0p0.dat")   .+ shiftPsi

    psiIM_ETL9p0ND   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-9p0.dat")   .+ shiftPsi .+ (9 * 0.04)/2
    psiIM_INTR9p0ND  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-9p0.dat")  .+ shiftPsi .+ (9 * 0.04)/2
    psiIM_HTL9p0ND   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-9p0.dat")   .+ shiftPsi .+ (9 * 0.04)/2

    psiIM_ETL24p0ND  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-24p0.dat")  .+ shiftPsi .+ (24 * 0.04)/2
    psiIM_INTR24p0ND = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-24p0.dat") .+ shiftPsi .+ (24 * 0.04)/2
    psiIM_HTL24p0ND  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-24p0.dat")  .+ shiftPsi .+ (24 * 0.04)/2

    psiIM_ETL30p0ND  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-30p0.dat")  .+ shiftPsi .+ (30 * 0.04)/2
    psiIM_INTR30p0ND = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-30p0.dat") .+ shiftPsi .+ (30 * 0.04)/2
    psiIM_HTL30p0ND  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-30p0.dat")  .+ shiftPsi .+ (30 * 0.04)/2
    ######################################################
    naIM_0p0ND       = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-0p0.dat")
    naIM_9p0ND       = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-9p0.dat")
    naIM_24p0ND      = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-24p0.dat")
    naIM_30p0ND      = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-30p0.dat")
    ######################################################
    dataND, sol0p0ND, sol9p0ND, sol24p0ND, sol30p0ND = PSCVolumeExclusion.main(test = true, plotting = false, comparison = true,
                                                                               eps = eps, nonlinearDiffusion = true)
    ipsi = dataND.index_psi
    ##############################################################################################
    ############################          modified drift data         ############################
    ##############################################################################################

    folder         = "drift"

    psiIM_ETL0p0MD   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-0p0.dat")   .+ shiftPsi
    psiIM_INTR0p0MD  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-0p0.dat")  .+ shiftPsi
    psiIM_HTL0p0MD   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-0p0.dat")   .+ shiftPsi

    psiIM_ETL9p0MD   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-9p0.dat")   .+ shiftPsi .+ (9 * 0.04)/2
    psiIM_INTR9p0MD  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-9p0.dat")  .+ shiftPsi .+ (9 * 0.04)/2
    psiIM_HTL9p0MD   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-9p0.dat")   .+ shiftPsi .+ (9 * 0.04)/2

    psiIM_ETL24p0MD  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-24p0.dat")  .+ shiftPsi .+ (24 * 0.04)/2
    psiIM_INTR24p0MD = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-24p0.dat") .+ shiftPsi .+ (24 * 0.04)/2
    psiIM_HTL24p0MD  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-24p0.dat")  .+ shiftPsi .+ (24 * 0.04)/2

    psiIM_ETL30p0MD  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-30p0.dat")  .+ shiftPsi .+ (30 * 0.04)/2
    psiIM_INTR30p0MD = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-30p0.dat") .+ shiftPsi .+ (30 * 0.04)/2
    psiIM_HTL30p0MD  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-30p0.dat")  .+ shiftPsi .+ (30 * 0.04)/2
    ######################################################
    naIM_0p0MD       = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-0p0.dat")
    naIM_9p0MD       = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-9p0.dat")
    naIM_24p0MD      = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-24p0.dat")
    naIM_30p0MD      = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-30p0.dat")
    ######################################################
    dataMD, sol0p0MD, sol9p0MD, sol24p0MD, sol30p0MD = PSCVolumeExclusion.main(test = true, plotting = false, comparison = true,
                                                                               eps = eps, nonlinearDiffusion = false)

    ##############################################################################################
    ############################         Calculate the errors         ############################
    ##############################################################################################

    ## w.r.t. electric potential
    psiErrCT0p0       = sol0p0MD[ipsi, :]  .- sol0p0ND[ipsi, :]
    psiErrCT9p0       = sol9p0MD[ipsi, :]  .- sol9p0ND[ipsi, :]
    psiErrCT24p0      = sol24p0MD[ipsi, :] .- sol24p0ND[ipsi, :]
    psiErrCT30p0      = sol30p0MD[ipsi, :] .- sol30p0ND[ipsi, :]

    psiErrIM_ETL0p0   = psiIM_ETL0p0MD     .- psiIM_ETL0p0ND
    psiErrIM_INTR0p0  = psiIM_INTR0p0MD    .- psiIM_INTR0p0ND
    psiErrIM_HTL0p0   = psiIM_HTL0p0MD     .- psiIM_HTL0p0ND

    psiErrIM_ETL9p0   = psiIM_ETL9p0MD     .- psiIM_ETL9p0ND
    psiErrIM_INTR9p0  = psiIM_INTR9p0MD    .- psiIM_INTR9p0ND
    psiErrIM_HTL9p0   = psiIM_HTL9p0MD     .- psiIM_HTL9p0ND

    psiErrIM_ETL24p0  = psiIM_ETL24p0MD    .- psiIM_ETL24p0ND
    psiErrIM_INTR24p0 = psiIM_INTR24p0MD   .- psiIM_INTR24p0ND
    psiErrIM_HTL24p0  = psiIM_HTL24p0MD    .- psiIM_HTL24p0ND

    psiErrIM_ETL30p0  = psiIM_ETL30p0MD    .- psiIM_ETL30p0ND
    psiErrIM_INTR30p0 = psiIM_INTR30p0MD   .- psiIM_INTR30p0ND
    psiErrIM_HTL30p0  = psiIM_HTL30p0MD    .- psiIM_HTL30p0ND

    ##############################################
    ##############################################
    na0p0ND     = get_density(sol0p0ND,  dataND, iphia, regionIntrinsic, inode = length_n+1:length_ni)
    na9p0ND     = get_density(sol9p0ND,  dataND, iphia, regionIntrinsic, inode = length_n+1:length_ni)
    na24p0ND    = get_density(sol24p0ND, dataND, iphia, regionIntrinsic, inode = length_n+1:length_ni)
    na30p0ND    = get_density(sol30p0ND, dataND, iphia, regionIntrinsic, inode = length_n+1:length_ni)

    na0p0MD     = get_density(sol0p0MD,  dataMD, iphia, regionIntrinsic, inode = length_n+1:length_ni)
    na9p0MD     = get_density(sol9p0MD,  dataMD, iphia, regionIntrinsic, inode = length_n+1:length_ni)
    na24p0MD    = get_density(sol24p0MD, dataMD, iphia, regionIntrinsic, inode = length_n+1:length_ni)
    na30p0MD    = get_density(sol30p0MD, dataMD, iphia, regionIntrinsic, inode = length_n+1:length_ni)

    naErrCT0p0  = na0p0MD     .- na0p0ND
    naErrCT9p0  = na9p0MD     .- na9p0ND
    naErrCT24p0 = na24p0MD    .- na24p0ND
    naErrCT30p0 = na30p0MD    .- na30p0ND

    naErrIM0p0  = naIM_0p0MD  .- naIM_0p0ND
    naErrIM9p0  = naIM_9p0MD  .- naIM_9p0ND
    naErrIM24p0 = naIM_24p0MD .- naIM_24p0ND
    naErrIM30p0 = naIM_30p0MD .- naIM_30p0ND

    ##############################################################################################
    ############################            Plot the errors           ############################
    ##############################################################################################

    ## Some plotting information
    function put_axvspan()
        PyPlot.axvspan(0.0,                        (h_ndoping)/μm,                       facecolor=[179/255 225/255 227/255], alpha = 0.6)
        PyPlot.axvspan(h_ndoping/μm,               (h_ndoping+h_intrinsic)/μm,           facecolor=[255/255 161/255 148/255], alpha = 0.6)
        PyPlot.axvspan((h_ndoping+h_intrinsic)/μm, (h_ndoping+h_intrinsic+h_pdoping)/μm, facecolor=[207/255 230/255 202/255], alpha = 0.6)
    end

    linewCT = 7
    linewIM = 5

    PyPlot.plot(coord./μm, psiErrCT0p0,  linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
    PyPlot.plot(coord./μm, psiErrCT9p0,  linewidth = linewCT, color = "blue",    label = "t = 9 s")
    PyPlot.plot(coord./μm, psiErrCT24p0, linewidth = linewCT, color = "royalblue",  label = "t = 24 s")
    PyPlot.plot(coord./μm, psiErrCT30p0, linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")
    ####################################################################
    PyPlot.plot(grid_ETL./μm,  psiErrIM_ETL0p0',   linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_intr./μm, psiErrIM_INTR0p0',  linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_HTL./μm,  psiErrIM_HTL0p0',   linewidth = linewIM, linestyle = ":", color = "black")
    ###############
    PyPlot.plot(grid_ETL./μm,  psiErrIM_ETL9p0',   linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_intr./μm, psiErrIM_INTR9p0',  linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_HTL./μm,  psiErrIM_HTL9p0',   linewidth = linewIM, linestyle = ":", color = "black")
    ###############
    PyPlot.plot(grid_ETL./μm,  psiErrIM_ETL24p0',  linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_intr./μm, psiErrIM_INTR24p0', linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_HTL./μm,  psiErrIM_HTL24p0',  linewidth = linewIM, linestyle = ":", color = "black")
    ###############
    PyPlot.plot(grid_ETL./μm,  psiErrIM_ETL30p0',  linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_intr./μm, psiErrIM_INTR30p0', linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_HTL./μm,  psiErrIM_HTL30p0',  linewidth = linewIM, linestyle = ":", color = "black")

    put_axvspan()
    PyPlot.grid()
    PyPlot.xlim(0.5e-1, 4.2e-1)
    if eps == 0.01
        PyPlot.ylim(-0.3e-3, 0.3e-3)
    elseif eps == 0.9
        PyPlot.ylim(-0.3e-1, 0.3e-1)
    end
    PyPlot.xlabel("space [\$ \\mu\\!\\!\$ m]", fontsize=18)
    PyPlot.ylabel(" \$ \\psi_{\\mathrm{MD}} - \\psi_{\\mathrm{ND}} \$ [V]", fontsize=18)
    PyPlot.tick_params(axis="both", labelsize=18)
    PyPlot.legend(fancybox = true, loc = "lower right", fontsize=15)
    PyPlot.tight_layout()

    if saveFig
        savefig("DD-parameter-eps-$textEps-psi-diff.pdf")
    end
    ##############################################
    ##############################################
    scale = Ca
    PyPlot.figure()
    PyPlot.plot(grid_intr./μm, naErrCT0p0./scale,  linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
    PyPlot.plot(grid_intr./μm, naErrCT9p0./scale,  linewidth = linewCT, color = "blue",    label = "t = 9 s")
    PyPlot.plot(grid_intr./μm, naErrCT24p0./scale, linewidth = linewCT, color = "royalblue",  label = "t = 24 s")
    PyPlot.plot(grid_intr./μm, naErrCT30p0./scale, linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")

    PyPlot.plot(grid_intr./μm, naErrIM0p0'./scale,  linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_intr./μm, naErrIM9p0'./scale,  linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_intr./μm, naErrIM24p0'./scale, linewidth = linewIM, linestyle = ":", color = "black")
    PyPlot.plot(grid_intr./μm, naErrIM30p0'./scale, linewidth = linewIM, linestyle = ":", color = "black")

    put_axvspan()
    PyPlot.grid()
    PyPlot.xlim(0.5e-1, 4.2e-1)
    if eps == 0.01
        PyPlot.ylim(-2.0e-3, 4.5e-3)
    elseif eps == 0.9
        PyPlot.ylim(-2.0e-1, 4.5e-1)
    end
    PyPlot.xlabel("space [\$ \\mu\\!\\!\$ m]", fontsize=18)
    PyPlot.ylabel(" \$ (n_{\\mathrm{a}, \\mathrm{MD}} - n_{\\mathrm{a}, \\mathrm{ND}})/C_\\mathrm{a} \$ [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=18)
    PyPlot.tick_params(axis="both", labelsize=18)
    PyPlot.legend(fancybox = true, loc = "upper center", fontsize=15)
    PyPlot.tight_layout()

    if saveFig
        savefig("DD-parameter-eps-$textEps-na-diff.pdf")
    end

end # main

end # module