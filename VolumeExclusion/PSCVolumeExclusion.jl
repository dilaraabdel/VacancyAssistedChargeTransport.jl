
#=
Influence of volume exclusion effects on charge transport in PSCs.
=#


module PSCVolumeExclusion

using ChargeTransport
using ExtendableGrids
using PyPlot
using DelimitedFiles
using LinearAlgebra

include("../parameters/Params_PSC_PCBM_MAPI_Pedot.jl")

# Eps set = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 0.95, 0.9999]

function main(;test = false, plotting = true, saveFig = false,
               eps = 0.01, nonlinearDiffusion = false,
               comparison = false, steadyStateCalc = false)

    PyPlot.rc("font", family="serif", size=15)
    PyPlot.rc("mathtext", fontset="dejavuserif")

    PyPlot.close("all")
    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    ## For all other parameters, check out the parameter template!!

    ## vacancy parameter
    Nanion           = Ca/eps
    Ea_i             = vacancy_energy(eps)

    ## contact voltage
    contactVoltage   = 1.2                  * V

    ## primary data for I-V scan protocol
    scanrate         = 0.04 * V/s
    number_tsteps    = 101
    endVoltage       = contactVoltage # bias goes until the given voltage
    tend             = endVoltage/scanrate
    tvalues          = range(0, stop = tend, length = number_tsteps)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    ## read in the grid from Ionmonger
    shift     = readdlm("data/IM-DD-parameter-grid-ETL.dat")[1]
    grid_ETL  = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-ETL.dat") .- shift)
    length_n  = length(grid_ETL)
    grid_intr = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-intr.dat").- shift)
    length_ni = length_n + length(grid_intr)
    grid_HTL  = 1.0e-9 .* (readdlm("data/IM-DD-parameter-grid-HTL.dat") .- shift)

    coord     = glue(grid_ETL[:, 1], grid_intr[:, 1])
    coord     = glue(coord[:, 1], grid_HTL[:, 1])
    grid      = simplexgrid(coord)

    ## cellmask! for defining the subregions and assigning region number
    cellmask!(grid, [0.0],                    [h_ndoping],               regionDonor,     tol = 1.0e-18) # n-doped region   = 1
    cellmask!(grid, [h_ndoping],              [h_ndoping+h_intrinsic],   regionIntrinsic, tol = 1.0e-18) # intrinsic region = 2
    cellmask!(grid, [h_ndoping+h_intrinsic],  [h_total],                 regionAcceptor,  tol = 1.0e-18) # p-doped region   = 3

    bfacemask!(grid, [h_ndoping],             [h_ndoping],               bregionJ1)  # first  inner interface
    bfacemask!(grid, [h_ndoping+h_intrinsic], [h_ndoping + h_intrinsic], bregionJ2)  # second inner interface

    if test == false
        println("*** done\n")
    end


    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    ## Initialize Data instance and fill in predefined data
    data                               = Data(grid, numberOfCarriers)
    data.modelType                     = Transient
    data.F                             = [Boltzmann, Boltzmann, FermiDiracMinusOne]

    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                 bulk_recomb_Auger = true,
                                                                 bulk_recomb_radiative = true,
                                                                 bulk_recomb_SRH = true)
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.boundaryType[bregionDonor]    = OhmicContact
    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

    if nonlinearDiffusion
        data.fluxApproximation        .= DiffusionEnhanced
    else
        data.fluxApproximation         = [DiffusionEnhanced, DiffusionEnhanced, DiffusionEnhancedModifiedDrift]
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                              = Params(grid, numberOfCarriers)

    params.temperature                                  = T
    params.UT                                           = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                         = zn
    params.chargeNumbers[iphip]                         = zp
    params.chargeNumbers[iphia]                         = za

    for ireg in 1:numberOfRegions # interior region data

        params.dielectricConstant[ireg]                 = ε[ireg] * ε0

        ## effective DOS, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nn[ireg]
        params.densityOfStates[iphip, ireg]             = Np[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = En[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = Ep[ireg]

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]
        params.mobility[iphia, ireg]                    = μa[ireg]

        ## recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, data, EI[ireg])

    end

    # vacancy parameter
    params.densityOfStates[iphia, regionIntrinsic]      = Nanion
    params.bandEdgeEnergy[ iphia, regionIntrinsic]      = Ea_i

    ## interior doping
    params.doping[iphin, regionDonor]                   = Cn
    params.doping[iphia, regionIntrinsic]               = Ca
    params.doping[iphip, regionAcceptor]                = Cp

    data.params                                         = params
    ctsys                                               = System(grid, data, unknown_storage=:sparse)
    ipsi                                                = data.index_psi
    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control              = SolverControl()
    control.verbose      = false#true
    control.maxiters     = 300
    control.abstol       = 1.0e-6
    control.reltol       = 1.0e-6
    control.tol_round    = 1.0e-5
    control.max_round    = 4
    control.damp_initial = 0.5
    control.damp_growth  = 1.21 # >= 1

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    ## initialize solution and starting vectors
    inival  = unknowns(ctsys); sol     = unknowns(ctsys)
    sol0p0  = unknowns(ctsys); sol3p0  = unknowns(ctsys)
    sol9p0  = unknowns(ctsys); sol15p0 = unknowns(ctsys)
    sol24p0 = unknowns(ctsys); sol30p0 = unknowns(ctsys)

    sol0p0  = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    inival  = sol0p0

    if test == false
        println("*** done\n")
    end

    ################################################################################
        if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    ## for saving I-V data
    IV = zeros(0); biasValues = zeros(0)

    for istep = 2:number_tsteps

        t  = tvalues[istep]       # Actual time
        Δu = t * scanrate         # Applied voltage
        Δt = t - tvalues[istep-1] # Time step size

        set_contact!(ctsys, bregionAcceptor, Δu = Δu)

        if test == false
            println("time value: t = $(t)")
        end

        sol = solve(ctsys; inival = inival, tstep = Δt, control = control)

        ## get I-V data
        current = get_current_val(ctsys, sol, inival, Δt)
        push!(IV, current); push!(biasValues, Δu)

        if t == 3
            sol3p0  = sol
        elseif t == 9
            sol9p0  = sol
        elseif t == 15
            sol15p0 = sol
        elseif t == 24
            sol24p0 = sol
        elseif t == 30
            sol30p0 = sol
        end

        inival = sol

    end # time loop

    if comparison
        return data, sol0p0, sol9p0, sol24p0, sol30p0
    end


    # if saveData
    #     res = [biasValues IV]
    #     writedlm("data/$textEps/$folder/jl-IV-forward-$folder-eps-$textEps-no-generation.dat", res)
    # end

    if test == false
        println("*** done\n")
    end

    if steadyStateCalc
        ################################################################################
        if test == false
            println("Calculate Steady State")
        end
        ################################################################################

        number_tstepsStead  = 301
        tvaluesStead        = range(tend, stop = 130, length = number_tstepsStead)

        ## initialize Array saving solutions for each time step
        solphin       = zeros(Float64, length(coord), length(tvaluesStead))
        solphip       = zeros(Float64, length(coord), length(tvaluesStead))
        solphia       = zeros(Float64, length(coord), length(tvaluesStead))
        solpsi        = zeros(Float64, length(coord), length(tvaluesStead))

        solphin[:, 1] = inival[iphin, :]; solphip[:, 1] = inival[iphip, :]
        solphia[:, 1] = inival[iphia, :]; solpsi[ :, 1] = inival[ipsi,  :]

        for istep = 2:number_tstepsStead

            t  = tvaluesStead[istep]       # Actual time
            Δt = t - tvaluesStead[istep-1] # Time step size

            if test == false
                println("          time value: Δt = $(t)")
            end

            sol    = solve(ctsys; inival = inival, tstep = Δt, control = control)
            inival = sol

            solphin[:, istep] .= sol[iphin, :]; solphip[:, istep] .= sol[iphip, :]
            solphia[:, istep] .= sol[iphia, :]; solpsi[ :, istep] .= sol[ipsi,  :]

        end

        # save l2 errors
        totall2Err   = zeros(0); psil2Err       = zeros(0)
        iphinl2Err   = zeros(0); iphipl2Err     = zeros(0)
        iphial2Err   = zeros(0)

        solStead     = [solphin[:, end] solphip[:, end] solphia[:, end] solpsi[:, end]]'

        for istep = 1:length(tvaluesStead)

            sol = [solphin[:, istep] solphip[:, istep] solphia[:, istep] solpsi[:, istep]]'

            err_psi = 0.0; err_n = 0.0; err_p = 0.0; err_a = 0.0 # L2 error

            for ix = 1:length(coord)-1
                h       = coord[ix+1] - coord[ix]
                # evaluation of L2 error contributions
                err_psi = err_psi + h * abs(sol[ipsi,  ix] - solStead[ipsi,  ix])^2
                err_n   = err_n   + h * abs(sol[iphin, ix] - solStead[iphin, ix])^2
                err_p   = err_p   + h * abs(sol[iphip, ix] - solStead[iphip, ix])^2
                if ix in length_n+1:length_ni-1
                    err_a   = err_a   + h * abs(sol[iphia, ix] - solStead[iphia, ix])^2
                end
            end

            err_tot = err_psi + err_n + err_p + err_a

            push!(totall2Err, err_tot); push!(psil2Err,   err_psi); push!(iphinl2Err, err_n)
            push!(iphipl2Err, err_p);   push!(iphial2Err, err_a)

        end

        return tvaluesStead, psil2Err, iphinl2Err, iphipl2Err, iphial2Err

        if test == false
            println("*** done\n")
        end

    end


    if plotting
        ################################################################################
        if test == false
            println("Some plotting")
        end
        ################################################################################

        ## read in measurement curve
        IV_meas = readdlm("data/IV-measured-CaladoNature2016-forward-dark.dat")

        ## Some plotting information
        function put_axvspan()
            PyPlot.axvspan(0.0,                        (h_ndoping)/μm,                       facecolor=[179/255 225/255 227/255], alpha = 0.6)
            PyPlot.axvspan(h_ndoping/μm,               (h_ndoping+h_intrinsic)/μm,           facecolor=[255/255 161/255 148/255], alpha = 0.6)
            PyPlot.axvspan((h_ndoping+h_intrinsic)/μm, (h_ndoping+h_intrinsic+h_pdoping)/μm, facecolor=[207/255 230/255 202/255], alpha = 0.6)
        end

        linewCT = 7; linewIM = 5

        #####################################################################################################
        #####################################################################################################

        plotIM = check_IMData_availability(eps) # check if Ionmonger data is available

        ## read in Ionmonger data
        if plotIM

            helpEps      = collect(string(eps)); helpEps[2] = 'p'
            textEps      = join(helpEps)

            if nonlinearDiffusion
                folder   = "diffusion"
            else
                folder   = "drift"
            end

            ## built-in voltage
            UT           = kB * T / q
            Vbi          = (En[regionDonor] - Ep[regionAcceptor] - kB * T * log(Nn[regionDonor] * Np[regionAcceptor]/ (Cn * Cp)))/q
            Nintr_d      = sqrt(Nn[regionDonor] * Np[regionDonor]) * exp(-(En[regionDonor] - Ep[regionDonor] ) / (2 * kB * T))
            psi0_d       = (En[regionDonor]+ Ep[regionDonor])/(2 * q) - 0.5 * UT * log(Nn[regionDonor]/Np[regionDonor]) + UT * asinh(Cn/(2*Nintr_d))
            shiftPsi     = - Vbi/2 + psi0_d

            if eps >= 0.999

                psi_im   = readdlm("data/no-ions/CT-DD-parameter-psi-t-0p0.dat")
                psi_im9  = readdlm("data/no-ions/CT-DD-parameter-psi-t-9p0.dat")
                psi_im24 = readdlm("data/no-ions/CT-DD-parameter-psi-t-24p0.dat")
                psi_im30 = readdlm("data/no-ions/CT-DD-parameter-psi-t-30p0.dat")

                ######################################################

                n_im         = readdlm("data/no-ions/CT-DD-parameter-n-intr-t-0p0.dat")'
                p_im         = readdlm("data/no-ions/CT-DD-parameter-p-intr-t-0p0.dat")'

                n_im9        = readdlm("data/no-ions/CT-DD-parameter-n-intr-t-9p0.dat")'
                p_im9        = readdlm("data/no-ions/CT-DD-parameter-p-intr-t-9p0.dat")'

                n_im24       = readdlm("data/no-ions/CT-DD-parameter-n-intr-t-24p0.dat")'
                p_im24       = readdlm("data/no-ions/CT-DD-parameter-p-intr-t-24p0.dat")'

                n_im30       = readdlm("data/no-ions/CT-DD-parameter-n-intr-t-30p0.dat")'
                p_im30       = readdlm("data/no-ions/CT-DD-parameter-p-intr-t-30p0.dat")'

                IV_im        = readdlm("data/no-ions/CT-DD-parameter-J.dat")
            else

                psiETL_im    = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-0p0.dat")   .+ shiftPsi
                psiINTR_im   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-0p0.dat")  .+ shiftPsi
                psiHTL_im    = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-0p0.dat")   .+ shiftPsi

                psiETL_im9   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-9p0.dat")   .+ shiftPsi .+ (9 * 0.04)/2
                psiINTR_im9  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-9p0.dat")  .+ shiftPsi .+ (9 * 0.04)/2
                psiHTL_im9   = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-9p0.dat")   .+ shiftPsi .+ (9 * 0.04)/2

                psiETL_im24  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-24p0.dat")  .+ shiftPsi .+ (24 * 0.04)/2
                psiINTR_im24 = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-24p0.dat") .+ shiftPsi .+ (24 * 0.04)/2
                psiHTL_im24  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-24p0.dat")  .+ shiftPsi .+ (24 * 0.04)/2

                psiETL_im30  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-ETL-t-30p0.dat")  .+ shiftPsi .+ (30 * 0.04)/2
                psiINTR_im30 = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-intr-t-30p0.dat") .+ shiftPsi .+ (30 * 0.04)/2
                psiHTL_im30  = readdlm("data/$textEps/$folder/IM-DD-parameter-psi-HTL-t-30p0.dat")  .+ shiftPsi .+ (30 * 0.04)/2

                ######################################################

                n_im         = readdlm("data/$textEps/$folder/IM-DD-parameter-n-intr-t-0p0.dat")
                p_im         = readdlm("data/$textEps/$folder/IM-DD-parameter-p-intr-t-0p0.dat")
                a_im         = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-0p0.dat")

                n_im9        = readdlm("data/$textEps/$folder/IM-DD-parameter-n-intr-t-9p0.dat")
                p_im9        = readdlm("data/$textEps/$folder/IM-DD-parameter-p-intr-t-9p0.dat")
                a_im9        = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-9p0.dat")

                n_im24       = readdlm("data/$textEps/$folder/IM-DD-parameter-n-intr-t-24p0.dat")
                p_im24       = readdlm("data/$textEps/$folder/IM-DD-parameter-p-intr-t-24p0.dat")
                a_im24       = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-24p0.dat")

                n_im30       = readdlm("data/$textEps/$folder/IM-DD-parameter-n-intr-t-30p0.dat")
                p_im30       = readdlm("data/$textEps/$folder/IM-DD-parameter-p-intr-t-30p0.dat")
                a_im30       = readdlm("data/$textEps/$folder/IM-DD-parameter-a-intr-t-30p0.dat")

                IV_im        = readdlm("data/$textEps/$folder/IM-DD-parameter-J.dat")
            end

        end

    end

    if plotting
        PyPlot.plot(coord./μm, sol0p0[ipsi, :],   linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
        PyPlot.plot(coord./μm, sol9p0[ipsi, :],   linewidth = linewCT, color = "blue",    label = "t = 9 s")
        PyPlot.plot(coord./μm, sol24p0[ipsi, :],  linewidth = linewCT, color = "royalblue",  label = "t = 24 s")
        PyPlot.plot(coord./μm, sol30p0[ipsi, :],  linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")
        ####################################################################
        if plotIM && eps < 0.999
            PyPlot.plot(grid_ETL./μm,  psiETL_im',    linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(grid_intr./μm, psiINTR_im',   linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(grid_HTL./μm,  psiHTL_im',    linestyle = ":", linewidth = linewIM, color = "black")
            ###############
            PyPlot.plot(grid_ETL./μm,  psiETL_im9',   linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(grid_intr./μm, psiINTR_im9',  linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(grid_HTL./μm,  psiHTL_im9',   linestyle = ":", linewidth = linewIM, color = "black")
            ###############
            PyPlot.plot(grid_ETL./μm,  psiETL_im24',  linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(grid_intr./μm, psiINTR_im24', linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(grid_HTL./μm,  psiHTL_im24',  linestyle = ":", linewidth = linewIM, color = "black")
            ###############
            PyPlot.plot(grid_ETL./μm,  psiETL_im30',  linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(grid_intr./μm, psiINTR_im30', linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(grid_HTL./μm,  psiHTL_im30',  linestyle = ":", linewidth = linewIM, color = "black")
        else
            PyPlot.plot(coord./μm, psi_im,   linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(coord./μm, psi_im9,  linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(coord./μm, psi_im24, linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.plot(coord./μm, psi_im30, linestyle = ":", linewidth = linewIM, color = "black")
        end

        put_axvspan()
        PyPlot.grid()
        PyPlot.xlim(0.7e-1, 4.0e-1)
        PyPlot.ylim(-5.05, -3.76)
        PyPlot.xlabel("space [\$ \\mu\\!\\!\$ m]", fontsize=18)
        PyPlot.ylabel("potential [V]", fontsize=18)
        PyPlot.tick_params(axis="both", labelsize=18)
        PyPlot.legend(fancybox = true, loc = "lower center")
        PyPlot.tight_layout()

        if saveFig
            savefig("DD-parameter-$folder-eps-$textEps-psi.pdf")
        end
    end

    ########################################################################################
    ########################################################################################

    if plotting && eps < 0.999
        a   = get_density(sol0p0,  data, iphia, regionIntrinsic, inode = length_n+1:length_ni)
        a9  = get_density(sol9p0,  data, iphia, regionIntrinsic, inode = length_n+1:length_ni)
        a24 = get_density(sol24p0, data, iphia, regionIntrinsic, inode = length_n+1:length_ni)
        a30 = get_density(sol30p0, data, iphia, regionIntrinsic, inode = length_n+1:length_ni)


        PyPlot.figure()
        subplot(211)
        PyPlot.semilogy(grid_intr./μm, a,   linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
        PyPlot.semilogy(grid_intr./μm, a9,  linewidth = linewCT, color = "blue",         label = "t = 9 s")
        PyPlot.semilogy(grid_intr./μm, a24, linewidth = linewCT, color = "royalblue",    label = "t = 24 s")
        PyPlot.semilogy(grid_intr./μm, a30, linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")
        ####################################################################
        if plotIM
            PyPlot.semilogy(grid_intr./μm, a_im',   linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, a_im9',  linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, a_im24', linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, a_im30', linestyle = ":", linewidth = linewIM, color = "black")
        end

        put_axvspan()
        PyPlot.grid()
        PyPlot.xlim(0.80e-1, 1.2e-1)
        PyPlot.ylim(1.0e21, 1.0e26)
        PyPlot.yticks( [1.0e22, 1.0e24, 1.0e26] )
        PyPlot.xlabel("space [\$ \\mu \\!\\! \$ m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fancybox = true, loc = "lower right")
        PyPlot.tight_layout()

        subplot(212)
        PyPlot.semilogy(grid_intr./μm, a,   linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
        PyPlot.semilogy(grid_intr./μm, a9,  linewidth = linewCT, color = "blue",         label = "t = 9 s")
        PyPlot.semilogy(grid_intr./μm, a24, linewidth = linewCT, color = "royalblue",    label = "t = 24 s")
        PyPlot.semilogy(grid_intr./μm, a30, linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")
        ####################################################################
        if plotIM
            PyPlot.semilogy(grid_intr./μm, a_im',   linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, a_im9',  linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, a_im24', linestyle = ":", linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, a_im30', linestyle = ":", linewidth = linewIM, color = "black")
        end

        put_axvspan()
        PyPlot.grid()
        PyPlot.xlim(3.0e-1, 3.95e-1)
        PyPlot.ylim(6.0e22, 4.0e25)
        PyPlot.yticks( [6.0e22, 1.0e24, 4.0e25] )
        PyPlot.xlabel("space [\$ \\mu \\!\\! \$ m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fancybox = true, loc = "lower left")
        PyPlot.tight_layout()

        if saveFig
            savefig("DD-parameter-$folder-eps-$textEps-a.pdf")
        end
    end

    ########################################################################################
    ########################################################################################

    if plotting
        n   = get_density(sol0p0,  data, iphin, regionIntrinsic, inode = length_n+1:length_ni)
        n9  = get_density(sol9p0,  data, iphin, regionIntrinsic, inode = length_n+1:length_ni)
        n24 = get_density(sol24p0, data, iphin, regionIntrinsic, inode = length_n+1:length_ni)
        n30 = get_density(sol30p0, data, iphin, regionIntrinsic, inode = length_n+1:length_ni)

        PyPlot.figure()

        PyPlot.semilogy(grid_intr./μm, n,   linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
        PyPlot.semilogy(grid_intr./μm, n9,  linewidth = linewCT, color = "royalblue",    label = "t = 9 s")
        PyPlot.semilogy(grid_intr./μm, n24, linewidth = linewCT, color = "deepskyblue",  label = "t = 24 s")
        PyPlot.semilogy(grid_intr./μm, n30, linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")
        ####################################################################
        if plotIM
            PyPlot.semilogy(grid_intr./μm, n_im',   linestyle = ":",  linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, n_im9',  linestyle = ":",  linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, n_im24', linestyle = ":",  linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, n_im30', linestyle = ":",  linewidth = linewIM, color = "black")
        end

        put_axvspan()
        PyPlot.grid()
        PyPlot.xlim(0.7e-1, 4.0e-1)
        PyPlot.xlabel("space [m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fancybox = true, loc = "lower center")
        PyPlot.title("Electrons")
        PyPlot.tight_layout()
    end

    ########################################################################################
    ########################################################################################

    if plotting
        p   = get_density(sol0p0,  data, iphip, regionIntrinsic, inode = length_n+1:length_ni)
        p9  = get_density(sol9p0,  data, iphip, regionIntrinsic, inode = length_n+1:length_ni)
        p24 = get_density(sol24p0,  data, iphip, regionIntrinsic, inode = length_n+1:length_ni)
        p30 = get_density(sol30p0,  data, iphip, regionIntrinsic, inode = length_n+1:length_ni)

        PyPlot.figure()

        PyPlot.semilogy(grid_intr./μm, p,   linewidth = linewCT, color = "midnightblue", label = "t = 0 s")
        PyPlot.semilogy(grid_intr./μm, p9,  linewidth = linewCT, color = "royalblue",    label = "t = 9 s")
        PyPlot.semilogy(grid_intr./μm, p24, linewidth = linewCT, color = "deepskyblue",  label = "t = 24 s")
        PyPlot.semilogy(grid_intr./μm, p30, linewidth = linewCT, color = "lightskyblue", label = "t = 30 s")
        ####################################################################
        if plotIM
            PyPlot.semilogy(grid_intr./μm, p_im',   linestyle = ":",  linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, p_im9',  linestyle = ":",  linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, p_im24', linestyle = ":",  linewidth = linewIM, color = "black")
            PyPlot.semilogy(grid_intr./μm, p_im30', linestyle = ":", linewidth = linewIM, color = "black")
        end

        put_axvspan()
        PyPlot.grid()
        PyPlot.xlim(0.7e-1, 4.0e-1)
        PyPlot.xlabel("space [m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fancybox = true, loc = "lower center")
        PyPlot.title("Holes")
        PyPlot.tight_layout()
    end

    ########################################################################################
    ########################################################################################

    if plotting
        PyPlot.figure()
        PyPlot.plot(IV_meas[:, 1], IV_meas[:, 2],     linewidth = 3, color = "black",    label= "measured")
        PyPlot.plot(biasValues,    IV.*(cm^2).*1.0e3, linewidth = 3, color = "blue", label= "ChargeTransport", marker = "o", )
        if plotIM && eps < 0.999
            PyPlot.plot(biasValues,   -IV_im[2:end],      linewidth = 3, color = "red",    label= "Ionmonger",       linestyle = "--")
        elseif plotIM && eps >= 0.999
            PyPlot.plot(biasValues,   IV_im.*(cm^2).*1.0e3,      linewidth = 3, color = "red",    label= "Ionmonger",       linestyle = "--")
        end
        PyPlot.grid()
        PyPlot.xlabel("bias [V]")
        PyPlot.ylabel("current density [mA \$ \\mathrm{cm}^{-2}\$ ]")
        PyPlot.legend()
        PyPlot.xlim(-0.09, 1.25)
        PyPlot.ylim(-0.2, 7.0)
        PyPlot.tight_layout()

        if test == false
            println("*** done\n")
        end
    end



    testval = sum(filter(!isnan, sol30p0))/length(sol30p0) # when using sparse storage, we get NaN values in solution
    return testval



end #  main


function test()
    testval = -0.601863712139032
  main(test = true) ≈ testval
end

end # module

