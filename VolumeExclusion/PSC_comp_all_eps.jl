module PSC_comp_all_eps

using PyPlot
using DelimitedFiles
using ChargeTransport
using ExtendableGrids

include("../parameters/Params_PSC_PCBM_MAPI_Pedot.jl")
include("PSCVolumeExclusion.jl")
include("PSC_comp_one_eps.jl")

function main(;saveFig = false)

    PyPlot.rc("font", family="serif", size=15)
    PyPlot.rc("mathtext", fontset="dejavuserif")

    ## read in the grid from Ionmonger
    shift     = readdlm("data/IM-DD-parameter-grid-ETL.dat")[1]
    grid_ETL  = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-ETL.dat") .- shift)
    length_n  = length(grid_ETL)
    grid_intr = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-intr.dat").- shift)
    length_ni = length_n + length(grid_intr)
    grid_HTL  = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-HTL.dat") .- shift)

    coord     = glue(grid_ETL[:, 1], grid_intr[:, 1])
    coord     = glue(coord[:, 1], grid_HTL[:, 1])

    epsVec       = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9]
    PsiErrCT0p0  = zeros(0); PsiErrIM0p0  = zeros(0)
    PsiErrCT9p0  = zeros(0); PsiErrIM9p0  = zeros(0)
    PsiErrCT24p0 = zeros(0); PsiErrIM24p0 = zeros(0)
    PsiErrCT30p0 = zeros(0); PsiErrIM30p0 = zeros(0)
    ########################################
    naErrCT0p0   = zeros(0); naErrIM0p0   = zeros(0)
    naErrCT9p0   = zeros(0); naErrIM9p0   = zeros(0)
    naErrCT24p0  = zeros(0); naErrIM24p0  = zeros(0)
    naErrCT30p0  = zeros(0); naErrIM30p0  = zeros(0)

    ## built-in voltage
    UT           = kB * T / q
    Vbi          = (En[regionDonor] - Ep[regionAcceptor] - kB * T * log(Nn[regionDonor] * Np[regionAcceptor]/ (Cn * Cp)))/q
    Nintr_d      = sqrt(Nn[regionDonor] * Np[regionDonor]) * exp(-(En[regionDonor] - Ep[regionDonor] ) / (2 * kB * T))
    psi0_d       = (En[regionDonor]+ Ep[regionDonor])/(2 * q) - 0.5 * UT * log(Nn[regionDonor]/Np[regionDonor]) + UT * asinh(Cn/(2*Nintr_d))
    shiftPsi     = - Vbi/2 + psi0_d

    ## loop over all eps values
    for eps in epsVec

        println("     eps = $eps")

        helpEps          = collect(string(eps)); helpEps[2] = 'p'
        textEps          = join(helpEps)
        ###################################################################
        ###############      nonlinear diffusion data       ###############
        ###################################################################

        folder           = "diffusion"

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
        ###################################################################
        ###############         modified drift data         ###############
        ###################################################################
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

        ###################################################################
        ###############           calculate error           ###############
        ###################################################################

        PsierrCT0p0       = findmax(abs.(sol0p0MD[ipsi, :]  .- sol0p0ND[ipsi, :]))[1]
        PsierrCT9p0       = findmax(abs.(sol9p0MD[ipsi, :]  .- sol9p0ND[ipsi, :]))[1]
        PsierrCT24p0      = findmax(abs.(sol24p0MD[ipsi, :] .- sol24p0ND[ipsi, :]))[1]
        PsierrCT30p0      = findmax(abs.(sol30p0MD[ipsi, :] .- sol30p0ND[ipsi, :]))[1]

        push!(PsiErrCT0p0, PsierrCT0p0); push!(PsiErrCT24p0, PsierrCT24p0)
        push!(PsiErrCT9p0, PsierrCT9p0); push!(PsiErrCT30p0, PsierrCT30p0)
        ######################################################
        PsiErrIM_ETL0p0   = findmax(abs.(psiIM_ETL0p0MD     .- psiIM_ETL0p0ND))[1]
        PsiErrIM_INTR0p0  = findmax(abs.(psiIM_INTR0p0MD    .- psiIM_INTR0p0ND))[1]
        PsiErrIM_HTL0p0   = findmax(abs.(psiIM_HTL0p0MD     .- psiIM_HTL0p0ND))[1]
        PsierrIM0p0       = findmax( [PsiErrIM_ETL0p0 PsiErrIM_INTR0p0 PsiErrIM_HTL0p0] )[1]
        ######################################################
        PsiErrIM_ETL9p0   = findmax(abs.(psiIM_ETL9p0MD     .- psiIM_ETL9p0ND))[1]
        PsiErrIM_INTR9p0  = findmax(abs.(psiIM_INTR9p0MD    .- psiIM_INTR9p0ND))[1]
        PsiErrIM_HTL9p0   = findmax(abs.(psiIM_HTL9p0MD     .- psiIM_HTL9p0ND))[1]
        PsierrIM9p0       = findmax( [PsiErrIM_ETL9p0 PsiErrIM_INTR9p0 PsiErrIM_HTL9p0] )[1]
        ######################################################
        PsiErrIM_ETL24p0  = findmax(abs.(psiIM_ETL24p0MD    .- psiIM_ETL24p0ND))[1]
        PsiErrIM_INTR24p0 = findmax(abs.(psiIM_INTR24p0MD   .- psiIM_INTR24p0ND))[1]
        PsiErrIM_HTL24p0  = findmax(abs.(psiIM_HTL24p0MD    .- psiIM_HTL24p0ND))[1]
        PsierrIM24p0      = findmax( [PsiErrIM_ETL24p0 PsiErrIM_INTR24p0 PsiErrIM_HTL24p0] )[1]
        ######################################################
        PsiErrIM_ETL30p0  = findmax(abs.(psiIM_ETL30p0MD    .- psiIM_ETL30p0ND))[1]
        PsiErrIM_INTR30p0 = findmax(abs.(psiIM_INTR30p0MD   .- psiIM_INTR30p0ND))[1]
        PsiErrIM_HTL30p0  = findmax(abs.(psiIM_HTL30p0MD    .- psiIM_HTL30p0ND))[1]
        PsierrIM30p0      = findmax( [PsiErrIM_ETL30p0 PsiErrIM_INTR30p0 PsiErrIM_HTL30p0] )[1]

        push!(PsiErrIM0p0, PsierrIM0p0); push!(PsiErrIM24p0, PsierrIM24p0)
        push!(PsiErrIM9p0, PsierrIM9p0); push!(PsiErrIM30p0, PsierrIM30p0)
        ######################################################
        ######################################################
        na0p0ND     = get_density(sol0p0ND,  dataND, iphia, regionIntrinsic, inode = length_n+1:length_ni-1)
        na9p0ND     = get_density(sol9p0ND,  dataND, iphia, regionIntrinsic, inode = length_n+1:length_ni-1)
        na24p0ND    = get_density(sol24p0ND, dataND, iphia, regionIntrinsic, inode = length_n+1:length_ni-1)
        na30p0ND    = get_density(sol30p0ND, dataND, iphia, regionIntrinsic, inode = length_n+1:length_ni-1)

        na0p0MD     = get_density(sol0p0MD,  dataMD, iphia, regionIntrinsic, inode = length_n+1:length_ni-1)
        na9p0MD     = get_density(sol9p0MD,  dataMD, iphia, regionIntrinsic, inode = length_n+1:length_ni-1)
        na24p0MD    = get_density(sol24p0MD, dataMD, iphia, regionIntrinsic, inode = length_n+1:length_ni-1)
        na30p0MD    = get_density(sol30p0MD, dataMD, iphia, regionIntrinsic, inode = length_n+1:length_ni-1)

        naerrCT0p0  = findmax(abs.(na0p0MD     .- na0p0ND))[1]
        naerrCT9p0  = findmax(abs.(na9p0MD     .- na9p0ND))[1]
        naerrCT24p0 = findmax(abs.(na24p0MD    .- na24p0ND))[1]
        naerrCT30p0 = findmax(abs.(na30p0MD    .- na30p0ND))[1]

        push!(naErrCT0p0,  naerrCT0p0); push!(naErrCT24p0, naerrCT24p0)
        push!(naErrCT9p0,  naerrCT9p0); push!(naErrCT30p0, naerrCT30p0)
        ######################################################
        naerrIM0p0  = findmax(abs.(naIM_0p0MD  .- naIM_0p0ND))[1]
        naerrIM9p0  = findmax(abs.(naIM_9p0MD  .- naIM_9p0ND))[1]
        naerrIM24p0 = findmax(abs.(naIM_24p0MD .- naIM_24p0ND))[1]
        naerrIM30p0 = findmax(abs.(naIM_30p0MD .- naIM_30p0ND))[1]

        push!(naErrIM0p0,  naerrIM0p0); push!(naErrIM24p0, naerrIM24p0)
        push!(naErrIM9p0,  naerrIM9p0); push!(naErrIM30p0, naerrIM30p0)


    end # for loop

    println("loop done ****")
    ###################################################################
    ###############               Summary               ###############
    ###################################################################

    linewCT = 7
    linewIM = 5

    PyPlot.plot(epsVec, PsiErrCT0p0,  markersize = 15, marker = "o", linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
    PyPlot.plot(epsVec, PsiErrCT9p0,  markersize = 15, marker = "o", linewidth = linewCT, color = "blue",    label = "t = 9 s")
    PyPlot.plot(epsVec, PsiErrCT24p0, markersize = 15, marker = "o", linewidth = linewCT, color = "royalblue",  label = "t = 24 s")
    PyPlot.plot(epsVec, PsiErrCT30p0, markersize = 15, marker = "o", linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")
    ##################################################
    PyPlot.plot(epsVec, PsiErrIM0p0,  linestyle = ":", markersize = 10, marker = "x",linewidth = linewIM, color = "black")
    PyPlot.plot(epsVec, PsiErrIM9p0,  linestyle = ":", markersize = 10, marker = "x",linewidth = linewIM, color = "black")
    PyPlot.plot(epsVec, PsiErrIM24p0, linestyle = ":", markersize = 10, marker = "x",linewidth = linewIM, color = "black")
    PyPlot.plot(epsVec, PsiErrIM30p0, linestyle = ":", markersize = 10, marker = "x",linewidth = linewIM, color = "black")

    PyPlot.grid()
    PyPlot.xlabel(" \$ \\epsilon \$ ", fontsize=17)
    PyPlot.ylabel(" \$|| \\psi_{\\mathrm{MD}} - \\psi_{\\mathrm{ND}}  ||_{L^{\\infty}} \$", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize=16)
    PyPlot.tight_layout()

    if saveFig
        savefig("DD-parameter-psi-diff.pdf")
    end

    ######################################################
    ######################################################
    PyPlot.figure()
    scale = Ca
    PyPlot.plot(epsVec, naErrCT0p0./scale,  markersize = 15, marker = "o", linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
    PyPlot.plot(epsVec, naErrCT9p0./scale,  markersize = 15, marker = "o", linewidth = linewCT, color = "blue",    label = "t = 9 s")
    PyPlot.plot(epsVec, naErrCT24p0./scale, markersize = 15, marker = "o", linewidth = linewCT, color = "royalblue",  label = "t = 24 s")
    PyPlot.plot(epsVec, naErrCT30p0./scale, markersize = 15, marker = "o", linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")
    ##################################################
    PyPlot.plot(epsVec, naErrIM0p0./scale,  linestyle = ":", markersize = 10, marker = "x",linewidth = linewIM, color = "black")
    PyPlot.plot(epsVec, naErrIM9p0./scale,  linestyle = ":", markersize = 10, marker = "x",linewidth = linewIM, color = "black")
    PyPlot.plot(epsVec, naErrIM24p0./scale, linestyle = ":", markersize = 10, marker = "x",linewidth = linewIM, color = "black")
    PyPlot.plot(epsVec, naErrIM30p0./scale, linestyle = ":", markersize = 10, marker = "x",linewidth = linewIM, color = "black")

    PyPlot.grid()
    PyPlot.xlabel(" \$ \\epsilon \$ ", fontsize=17)
    PyPlot.ylabel(" \$||(n_{\\mathrm{a}, \\mathrm{MD}} - n_{\\mathrm{a}, \\mathrm{ND}})/C_\\mathrm{a} ||_{L^{\\infty}} \$", fontsize=18)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.legend(fontsize=16)
    PyPlot.tight_layout()

    if saveFig
        savefig("DD-parameter-a-diff.pdf")
    end

end # main

end # module