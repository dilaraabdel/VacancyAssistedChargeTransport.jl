

#=
# Large time behavior of a three-layer perovskite model ETL|Mapi|HTL with physically
realistic parameters, present photogeneration and time-independent boundary conditions.

=#

module PSC

using ChargeTransport
using ExtendableGrids
using PyPlot
using LinearAlgebra

include("../parameters/Params_PSC_TiO2_MAPI_Pedot.jl")

function main(;plotting = false, saveFig = false, test = false)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # For all parameters, we refer to the parameters template, which is included at the
    # beginning of this script.

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    coordn   = collect(0.0:hn:h_ndoping)
    lengthc1 = length(coordn)
    coord_i  = collect(h_ndoping:hi:h_ndoping+h_intrinsic)
    lengthc2 = lengthc1 + length(coord_i) - 1
    coord_p  = collect(h_ndoping+h_intrinsic:hp:h_total)

    coord    = glue(coordn, coord_i)
    coord    = glue(coord,   coord_p)
    grid     = simplexgrid(coord)

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

    # initialize Data instance and fill in data
    data                               = Data(grid, numberOfCarriers)
    data.modelType                     = Transient
    data.F                             = [Boltzmann, Boltzmann, FermiDiracMinusOne]
    data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                bulk_recomb_Auger = false,
                                                                bulk_recomb_radiative = true,
                                                                bulk_recomb_SRH = true)

    data.boundaryType[bregionDonor]    = OhmicContact
    data.boundaryType[bregionAcceptor] = OhmicContact
    data.generationModel               = GenerationBeerLambert
    data.fluxApproximation            .= ExcessChemicalPotential

    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionIntrinsic])

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

        # effective dos, band edge energy and mobilities
        params.densityOfStates[iphin, ireg]             = Nn[ireg]
        params.densityOfStates[iphip, ireg]             = Np[ireg]
        params.densityOfStates[iphia, ireg]             = Na[ireg]

        params.bandEdgeEnergy[iphin, ireg]              = En[ireg]
        params.bandEdgeEnergy[iphip, ireg]              = Ep[ireg]
        params.bandEdgeEnergy[iphia, ireg]              = Ea[ireg]

        params.mobility[iphin, ireg]                    = μn[ireg]
        params.mobility[iphip, ireg]                    = μp[ireg]
        params.mobility[iphia, ireg]                    = μa[ireg]

        # recombination parameters
        params.recombinationRadiative[ireg]             = r0[ireg]
        params.recombinationSRHLifetime[iphin, ireg]    = τn[ireg]
        params.recombinationSRHLifetime[iphip, ireg]    = τp[ireg]
        params.recombinationSRHTrapDensity[iphin, ireg] = trap_density!(iphin, ireg, data, EI[ireg])
        params.recombinationSRHTrapDensity[iphip, ireg] = trap_density!(iphip, ireg, data, EI[ireg])

        ## generation parameters
        params.generationIncidentPhotonFlux[ireg]       = incidentPhotonFlux[ireg]
        params.generationAbsorption[ireg]               = absorption[ireg]
    end

    # parameter which passes the shift information in the Beer-Lambert generation
    params.generationPeak                               = h_ndoping

    # doping
    params.doping[iphin,  regionDonor]                  = Cn
    params.doping[iphia,  regionIntrinsic]              = Ca
    params.doping[iphip,  regionAcceptor]               = Cp

    data.params = params
    ctsys       = System(grid, data, unknown_storage=:sparse)
    ipsi        = data.index_psi

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define control parameters for Newton solver")
    end
    ################################################################################

    control              = SolverControl()
    control.maxiters     = 300
    control.abstol       = 1.0e-5
    control.reltol       = 1.0e-7
    control.tol_round    = 1.0e-5
    control.max_round    = 5
    control.damp_initial = 0.5
    control.damp_growth  = 1.21

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Compute solution in thermodynamic equilibrium")
    end
    ################################################################################

    # initialize solution and starting vectors
    inival  = unknowns(ctsys); initialCond  = unknowns(ctsys)
    sol     = unknowns(ctsys); solStat      = unknowns(ctsys)
    inival .= 0.0;             initialCond .= 0.0

    solEQ   = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    inival  = solEQ

    if test == false
        println("*** done\n")
    end

    ################################################################################
     if test == false
        println("Going til open circuit voltage and putting generation on")
    end
    ################################################################################

    data.calculationType = OutOfEquilibrium

    ## these values are needed for putting the generation slightly on
    I                    = collect(length(bias):-1:0.0)
    LAMBDA               = 10 .^ (-I)

    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 1.0e30
    ctsys.fvmsys.boundary_values[iphia, bregionJ2]  = phiaOC

    for ii in eachindex(bias)

        set_contact!(ctsys, bregionAcceptor, Δu = bias[ii])
        data.λ2   = LAMBDA[ii+1] ## turn slowly generation on

        if test == false
            println("increase bias: Δu = $(bias[ii])")
            println("increase generation with λ2 = $(data.λ2)")
        end

        sol    = solve(ctsys; inival = inival, control = control)
        inival = sol

    end

    initialCond = sol; inival = initialCond


    if test == false
        println("*** done\n")
    end
    ################################################################################
    if test == false
        println("Calculate steady state")
    end
    ################################################################################

    # adjust newton parameters
    control.abstol    = 1.0e-10
    control.reltol    = 1.0e-10
    control.tol_round = 1.0e-10

    # fix the boundary value of phia (Δψ corresponds to the value the electric potential would have)
    phiaConst = kB * T/q * ( log(Ca/Na_i) -  log(1 - Ca/Na_i) ) - Ea_i/q - Δψ

    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 1.0e30
    ctsys.fvmsys.boundary_values[iphia, bregionJ2]  = phiaConst

    # reference steady state
    ## We do not make use of the stationary solution in the following, it is just used for a consistency check at the end
    ## to see, if steady state and stationary solution coincide
    solStat       = solve(ctsys, inival = inival, control = control)

    ## initialize Array saving solutions for each time step
    solphin       = zeros(Float64, length(coord), length(tvalues))
    solphip       = zeros(Float64, length(coord), length(tvalues))
    solphia       = zeros(Float64, length(coord), length(tvalues))
    solpsi        = zeros(Float64, length(coord), length(tvalues))

    solphin[:, 1] = inival[iphin, :]; solphip[:, 1] = inival[iphip, :]
    solphia[:, 1] = inival[iphia, :]; solpsi[ :, 1] = inival[ipsi,  :]

    for istep = 2:numbertsteps

        t  = tvalues[istep]       # Actual time
        Δt = t - tvalues[istep-1] # Time step size

        if test == false
            println("          time value: Δt = $(t)")
        end

        sol    = solve(ctsys; inival = inival, tstep = Δt, control = control)
        inival = sol

        solphin[:, istep] .= sol[iphin, :]; solphip[:, istep] .= sol[iphip, :]
        solphia[:, istep] .= sol[iphia, :]; solpsi[ :, istep] .= sol[ipsi,  :]

    end

    ################################################################################
    if test == false
        println("Calculate error and relative entropy")
    end
    ################################################################################

    regionNodes = [2:lengthc1, lengthc1+1:lengthc2, lengthc2+1:length(coord)-1] # first node will be calculated seperately
    p           = data.params

    #######################
    # save relative entropy w.r.t. steady state
    relEntropy  = zeros(0); relEntropyPsi = zeros(0); relEntropyn = zeros(0)
    relEntropyp = zeros(0); relEntropya   = zeros(0)

    # save L2 errors
    totall2Err  = zeros(0); psil2Err      = zeros(0); iphinl2Err  = zeros(0)
    iphipl2Err  = zeros(0); iphial2Err    = zeros(0)

    solStead    = [solphin[:, end] solphip[:, end] solphia[:, end] solpsi[:, end]]'

    for istep = 1:length(tvalues)

        sol = [solphin[:, istep] solphip[:, istep] solphia[:, istep] solpsi[:, istep]]'

        E_psi  = 0.0; E_n  = 0.0; E_p  = 0.0; E_a  = 0.0 # entropy w.r.t. steady state

        #### relative entropy calculations ####
        for ireg in eachindex(regionNodes)

            regNode = regionNodes[ireg] # extract the nodes of selected region
            inode   = regNode[1]-1
            Nn      = data.params.densityOfStates[iphin, ireg]
            Np      = data.params.densityOfStates[iphip, ireg]

            # calculate first values
            n1      = get_density(sol, data, iphin, ireg, inode = inode)
            p1      = get_density(sol, data, iphip, ireg, inode = inode)
            nStead1 = get_density(solStead, data, iphin, ireg, inode = inode)
            pStead1 = get_density(solStead, data, iphip, ireg, inode = inode)

            if ireg == 2 # only in intrinsic layer
                Na      = data.params.densityOfStates[iphia, ireg]
                a1      = get_density(sol, data, iphia, ireg, inode = inode)
                aStead1 = get_density(solStead, data, iphia, ireg, inode = inode)
            end

            for ix in regNode

                h            = coord[ix+1] - coord[ix]
                gradpsiStead = (solStead[ipsi, ix+1] - solStead[ipsi, ix]) / h
                gradpsi      = (sol[ipsi, ix+1] - sol[ipsi, ix]) / h

                n2      = get_density(sol, data, iphin, ireg, inode = ix)
                p2      = get_density(sol, data, iphip, ireg, inode = ix)
                nStead2 = get_density(solStead, data, iphin, ireg, inode = ix)
                pStead2 = get_density(solStead, data, iphip, ireg, inode = ix)

                # evaluation of entropy contributions
                E_psi = E_psi + 0.5 * data.params.dielectricConstant[ireg] * h * ( gradpsi - gradpsiStead )^2
                E_n   = E_n   + 0.5 * h * ( HElectric(p, n1, nStead1, iphin, ireg) + HElectric(p, n2, nStead2, iphin, ireg) )
                E_p   = E_p   + 0.5 * h * ( HElectric(p, p1, pStead1, iphip, ireg) + HElectric(p, p2, pStead2, iphip, ireg) )

                if ireg == regionIntrinsic # only in intrinsic layer
                    a2      = get_density(sol, data, iphia, ireg, inode = ix)
                    aStead2 = get_density(solStead, data, iphia, ireg, inode = ix)

                    E_a        = E_a + 0.5 * h * ( HIonic(p, a1, aStead1, iphia, ireg) + HIonic(p, a2, aStead2, iphia, ireg) )

                    a1      = a2; aStead1 = aStead2

                else
                    E_a        = E_a   + 0.0
                end

                n1 = n2; nStead1 = nStead2
                p1 = p2; pStead1 = pStead2

            end

        end

        push!(relEntropyPsi, E_psi); push!(relEntropyn, E_n); push!(relEntropyp, E_p)
        push!(relEntropya,   E_a);   push!(relEntropy,  E_psi + E_n + E_p + E_a)

        #### L2 error calculation ####
        # we calculate the L2 error on each subdomain
        err_psi =   sqrt(hn) .* norm( sol[ipsi,  1:lengthc1-1]      - solStead[ipsi,  1:lengthc1-1],      2)
                  + sqrt(hi) .* norm( sol[ipsi,  lengthc1:lengthc2] - solStead[ipsi,  lengthc1:lengthc2], 2)
                  + sqrt(hp) .* norm( sol[ipsi,  lengthc2+1:end]    - solStead[ipsi,  lengthc2+1:end],    2)
        ####################
        err_n   =   sqrt(hn) .* norm( sol[iphin, 1:lengthc1-1]      - solStead[iphin, 1:lengthc1-1],      2)
                  + sqrt(hi) .* norm( sol[iphin, lengthc1:lengthc2] - solStead[iphin, lengthc1:lengthc2], 2)
                  + sqrt(hp) .* norm( sol[iphin, lengthc2+1:end]    - solStead[iphin, lengthc2+1:end],    2)
        ####################
        err_p   =   sqrt(hn) .* norm( sol[iphip, 1:lengthc1-1]      - solStead[iphip, 1:lengthc1-1],      2)
                  + sqrt(hi) .* norm( sol[iphip, lengthc1:lengthc2] - solStead[iphip, lengthc1:lengthc2], 2)
                  + sqrt(hp) .* norm( sol[iphip, lengthc2+1:end]    - solStead[iphip, lengthc2+1:end],    2)
        ####################
        err_a   =   sqrt(hi) .* norm( sol[iphia, lengthc1:lengthc2] - solStead[iphia, lengthc1:lengthc2], 2)

        err_psi = norm( sol[ipsi, :]  - solStead[ipsi, :],  2 )
        err_n   = norm( sol[iphin, :] - solStead[iphin, :], 2 )
        err_p   = norm( sol[iphip, :] - solStead[iphip, :], 2 )
        err_a   = norm( sol[iphia, lengthc1:lengthc2] - solStead[iphia, lengthc1:lengthc2],  2 )

        err_tot = err_psi + err_n + err_p + err_a

        push!(totall2Err, err_tot); push!(psil2Err,   err_psi); push!(iphinl2Err, err_n)
        push!(iphipl2Err, err_p);   push!(iphial2Err, err_a)

    end

    if test == false
        println("*** done\n")
    end

    if plotting
        ################################################################################
        if test == false
            println("Plotting")
        end
        ################################################################################
        PyPlot.figure() ## Plot initial condition densities

        nn1 = get_density(sol, data, iphin, regionDonor, inode = 1:lengthc1)
        nn2 = get_density(sol, data, iphin, regionIntrinsic, inode = lengthc1:lengthc2)
        nn3 = get_density(sol, data, iphin, regionAcceptor, inode = lengthc2:length(coord))

        np1 = get_density(sol, data, iphip, regionDonor, inode = 1:lengthc1)
        np2 = get_density(sol, data, iphip, regionIntrinsic, inode = lengthc1:lengthc2)
        np3 = get_density(sol, data, iphip, regionAcceptor, inode = lengthc2:length(coord))

        na  = get_density(sol, data, iphia, regionIntrinsic, inode = lengthc1:lengthc2)
        ###########################
        nn1IC = get_density(initialCond, data, iphin, regionDonor, inode = 1:lengthc1)
        nn2IC = get_density(initialCond, data, iphin, regionIntrinsic, inode = lengthc1:lengthc2)
        nn3IC = get_density(initialCond, data, iphin, regionAcceptor, inode = lengthc2:length(coord))

        np1IC = get_density(initialCond, data, iphip, regionDonor, inode = 1:lengthc1)
        np2IC = get_density(initialCond, data, iphip, regionIntrinsic, inode = lengthc1:lengthc2)
        np3IC = get_density(initialCond, data, iphip, regionAcceptor, inode = lengthc2:length(coord))

        naIC  = get_density(initialCond, data, iphia, regionIntrinsic, inode = lengthc1:lengthc2)
        ###########################

        subplot(311)
        PyPlot.semilogy(coord[1:lengthc1]./μm,        nn1, color= "green", linewidth = 3, label =  "\$ n_{\\mathrm{n}}^\\infty \$" )
        PyPlot.semilogy(coord[lengthc1:lengthc2]./μm, nn2, color= "green", linewidth = 3 )
        PyPlot.semilogy(coord[lengthc2:end]./μm,      nn3, color= "green", linewidth = 3 )
        ##########
        PyPlot.semilogy(coord[1:lengthc1]./μm,        nn1IC, linewidth = 3, linestyle= ":", color="black", label = "IC" )
        PyPlot.semilogy(coord[lengthc1:lengthc2]./μm, nn2IC, linewidth = 3, linestyle= ":", color="black" )
        PyPlot.semilogy(coord[lengthc2:end]./μm,      nn3IC, linewidth = 3, linestyle= ":", color="black" )
        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=16)
        PyPlot.ylim(1.0e-3, 1.0e27)
        PyPlot.yticks( [1.0e-3, 1.0e12, 1.0e27] )
        PyPlot.tight_layout()
        ########################
        subplot(312)
        PyPlot.semilogy(coord[1:lengthc1]./μm,        np1, color= "red",  linewidth = 3, label =  "\$ n_{\\mathrm{p}}^\\infty\$" )
        PyPlot.semilogy(coord[lengthc1:lengthc2]./μm, np2, color= "red",  linewidth = 3)
        PyPlot.semilogy(coord[lengthc2:end]./μm,      np3, color= "red",  linewidth = 3 )
        ##########
        PyPlot.semilogy(coord[1:lengthc1]./μm,        np1IC, linewidth = 3, linestyle= ":", color="black" )
        PyPlot.semilogy(coord[lengthc1:lengthc2]./μm, np2IC, linewidth = 3, linestyle= ":", color="black" )
        PyPlot.semilogy(coord[lengthc2:end]./μm,      np3IC, linewidth = 3, linestyle= ":", color="black" )
        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=16)
        PyPlot.ylim(1.0e-3, 1.0e27)
        PyPlot.yticks( [1.0e-3, 1.0e12, 1.0e27] )
        PyPlot.tight_layout()
        ##########
        subplot(313)
        PyPlot.semilogy(coord[lengthc1:lengthc2]./μm, na,  color= "gold", linewidth = 3, label =  "\$ n_{\\mathrm{a}}^\\infty\$" )
        PyPlot.semilogy(coord[lengthc1:lengthc2]./μm, naIC, linewidth = 3, linestyle= ":", color="black" )
        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=16)
        PyPlot.xlim(-0.03, 0.73)
        PyPlot.ylim(1.0e23, 1.0e26)
        PyPlot.yticks( [1.0e23, 1.0e25, 1.0e26] )
        PyPlot.tight_layout()

        PyPlot.tight_layout()

        if saveFig
            savefig("steady-state-test-case-2-$paramsname-dens.pdf")
        end

        ##################################
        PyPlot.figure() ## Plot initial condition and steady state potentials

        subplot(211)
        PyPlot.plot(coord./μm, sol[iphin, :], color= "green", linestyle = "-", linewidth = 3, label =  "\$ \\varphi_{\\mathrm{n}}^{\\infty}\$" )
        PyPlot.plot(coord./μm, sol[iphip, :], color= "red",   linewidth = 3, label =  "\$ \\varphi_{\\mathrm{p}}^{\\infty}\$" )
        PyPlot.plot(coord./μm, sol[iphia, :], color= "gold", linewidth = 3, label =  "\$ \\varphi_{\\mathrm{a}}^{\\infty}\$" )

        PyPlot.plot(coord./μm, initialCond[iphin, :], linewidth = 3, linestyle= ":", color="black", label= "IC" )
        PyPlot.plot(coord./μm, initialCond[iphip, :], linewidth = 3, linestyle= ":", color="black")
        PyPlot.plot(coord[lengthc1:lengthc2]./μm, initialCond[iphia, lengthc1:lengthc2], linewidth = 3, linestyle= ":", color="black" )

        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("potential [V]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=14)
        if paramsname == "IM"
            PyPlot.ylim(-0.3, 1.1)
        end
        PyPlot.tight_layout()

        ####################
        subplot(212)
        PyPlot.plot(coord./μm, sol[ipsi,:], label = "\$\\psi^{\\infty}\$", color="b", linewidth= 3)
        PyPlot.plot(coord./μm, initialCond[ipsi, :], linewidth = 3, linestyle= ":", color="black", label= "IC")

        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("potential [V]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        if paramsname == "IM"
            PyPlot.ylim(-4.5, -4.0)
        end
        PyPlot.legend(fontsize=16)
        PyPlot.tight_layout()
        if saveFig
            savefig("steady-state-test-case-2-$paramsname.pdf")
        end

    end

    if plotting

        ################################################################################
        #  Putting a tolerance
        ################################################################################

        function adjustWithTolerance!(xx)
            eps = 5.0e-18
            for i in eachindex(xx)
                if xx[i] <= eps
                    xx[i] = eps
                end
            end
            return xx
        end

        if paramsname == "IM"
            plotStep = 15
        else
            plotStep = 8
        end

        PyPlot.figure()

        PyPlot.semilogy(tvalues[1:2:end], adjustWithTolerance!(relEntropy[1:2:end]),  color = "darkgreen", marker = "x")
        PyPlot.xlabel("time [s]", fontsize=17)
        PyPlot.ylabel("relative entropy w.r.t. steady state", fontsize=17)
        PyPlot.grid()
        PyPlot.tight_layout()
        if paramsname == "IM"
            PyPlot.xlim(-15, 230)
            PyPlot.ylim(0.2e-18, 6.0e1)
        else
            PyPlot.xlim(-30.0, 1200)
            PyPlot.ylim(0.2e-18, 6.0e1)
        end
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()

        if saveFig
            savefig("rel-entropy-test-case-2-$paramsname.pdf")
        end

        ######################################################

        PyPlot.figure()
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!((psil2Err[1:plotStep:end]).^2),   color = "blue",  marker = "x", linewidth =4, label = "\$ \\psi  - \\psi^\\infty \$")
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!((iphinl2Err[1:plotStep:end]).^2), color = "green", marker = "x", linewidth =4,  label = "\$ \\varphi_{\\mathrm{n}} - \\varphi_{\\mathrm{n}}^\\infty \$")
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!((iphipl2Err[1:plotStep:end]).^2), color = "red",  linestyle="dashed",  marker = "o", label = "\$ \\varphi_{\\mathrm{p}} - \\varphi_{\\mathrm{p}}^\\infty \$")
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!((iphial2Err[1:plotStep:end]).^2), color = "gold", linestyle="dashed", marker = "o", label = "\$ \\varphi_{\\mathrm{a}} - \\varphi_{\\mathrm{a}}^\\infty \$")

        PyPlot.xlabel("time [s]", fontsize=17)
        PyPlot.ylabel("\$ || u - u^{\\infty} ||_{L^2}^2\$", fontsize=17)
        PyPlot.legend(fontsize=16)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.grid()
        PyPlot.tight_layout()
        if paramsname == "IM"
            PyPlot.xlim(-15, 230)
            PyPlot.ylim(0.2e-18, 6.0e1)
        end
        if saveFig
            savefig("l2-error-test-case-2-$paramsname.pdf")
        end

    end

    testval = sum(filter(!isnan, solStat))/length(solStat) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
      testval = -0.775820946182112 # IM
    main(test = true) ≈ testval
end


end # module
