

#=
# Large time behavior of a TMDC memristive device with physically realistic parameters
and and time-independent boundary conditions.

=#

module TMDC

using ChargeTransport
using ExtendableGrids
using PyPlot
using LinearAlgebra
using Roots
using Trapz


include("../parameters/Params_TMDC_S1.jl")

function main(;plotting = true, saveFig = false, test = false,
               naAvgSweep = false,  EaSweep = -4.32*eV)


    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Define the scan protocol")
    end
    ################################################################################

    biasValQuasiStatic = 1.6 * V

    ## Define scan protocol function
    function scanProtocol(t)

        t0 = t % period

        if    0.0 <= t0  && t0 <= period/4
            biasVal = 0.0 + scanrate * t0
        elseif  t0 >= period/4  && t0 <= 3*period/4
            biasVal = amplitude .- scanrate *(t0-period/4)
        elseif  t0 >= 3*period/4 && t0 <= period
            biasVal = - amplitude .+ scanrate * (t0-3*period/4)
        else
            biasVal = 0.0
        end

        return biasVal

    end

    g(t)            = scanProtocol(t) - biasValQuasiStatic
    tEndQuasiStatic = find_zero(g, 0.0) + period
    ## Define scan protocol function

    function scanProtocolQuasiStatic(t)

        if t <= tEndQuasiStatic
            biasVal = scanProtocol(t)
        else
            biasVal = biasValQuasiStatic
        end

        return biasVal

    end

    # Apply zero voltage on left boundary and the given scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, scanProtocolQuasiStatic]

    if naAvgSweep
        Ea2 = EaSweep
    else
        ## naAvg = 1.017728681081517e17  /m^3  with  Ea = -4.88  *eV
        ## naAvg = 1.1815615845736299e20 /m^3  with  Ea = -4.695 *eV
        ## naAvg = 6.423415588996168e23  /m^3  with  Ea = -4.32  *eV
        ## naAvg = 1.01999084051074e26   /m^3  with  Ea = -3.98  *eV
        ## naAvg = 1.0340477150552678e27 /m^3  with  Ea = -3.38  *eV
        Ea2 = Ea
    end

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

    # non-uniform grid
    coord1       = geomspace(0.0,       h_flake/2, 3e-4 * h_flake, 2e-2 * h_flake)
    coord2       = geomspace(h_flake/2, h_flake,   2e-2 * h_flake, 3e-4 * h_flake)
    coord        = glue(coord1, coord2)

    grid         = simplexgrid(coord)

    ## set region in grid
    cellmask!(grid, [0.0], [h_flake], regionflake, tol = 1.0e-18)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define System and fill in information about model")
    end
    ################################################################################

    # initialize Data instance and fill in data
    data                            = Data(grid, numberOfCarriers, contactVoltageFunction = contactVoltageFunction)
    data.modelType                  = Transient
    data.F                          = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]
    data.bulkRecombination          = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                             bulk_recomb_Auger = false,
                                                             bulk_recomb_radiative = false,
                                                             bulk_recomb_SRH = false)

    data.boundaryType[bregionLeft]  = MixedOhmicSchottkyContact
    data.boundaryType[bregionRight] = MixedOhmicSchottkyContact
    data.fluxApproximation         .= ExcessChemicalPotential

    enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionflake])

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                     = Params(grid, numberOfCarriers)

    params.temperature                         = T
    params.UT                                  = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                = zn
    params.chargeNumbers[iphip]                = zp
    params.chargeNumbers[iphia]                = za

    params.dielectricConstant[regionflake]     = εr * ε0

    ## effective DOS, band-edge energy and mobilities
    params.densityOfStates[iphin, regionflake] = Nn
    params.densityOfStates[iphip, regionflake] = Np
    params.densityOfStates[iphia, regionflake] = Na

    params.bandEdgeEnergy[iphin, regionflake]  = En
    params.bandEdgeEnergy[iphip, regionflake]  = Ep
    params.bandEdgeEnergy[iphia, regionflake]  = Ea2

    params.mobility[iphin, regionflake]        = μn
    params.mobility[iphip, regionflake]        = μp
    params.mobility[iphia, regionflake]        = μa

    ## doping
    params.doping[iphin, regionflake]          = Cn

    ## in case of Schottky contact
    params.SchottkyBarrier[bregionLeft]        = barrierLeft
    params.SchottkyBarrier[bregionRight]       = barrierRight

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
    control.damp_initial = damp_initial
    control.damp_growth  = damp_growth
    control.max_round    = max_round
    control.maxiters     = maxiters

    control.abstol       = abstol
    control.reltol       = reltol
    control.tol_round    = tol_round

    control.Δu_opt       = Δu_opt
    control.Δt           = Δt
    control.Δt_min       = Δt_min
    control.Δt_max       = Δt_max
    control.Δt_grow      = Δt_grow

    if test == false
        control.verbose  = "ed"
    else
        control.verbose = ""
    end


    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define initial condition")
    end
    ################################################################################

    # initialize solution and starting vectors
    inival      = unknowns(ctsys); inival      .= 0.0
    initialCond = unknowns(ctsys); initialCond .= 0.0
    sol         = unknowns(ctsys);

    solEQ       = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Calculate steady state")
    end
    ################################################################################

    data.calculationType = OutOfEquilibrium

    # we go until the correct bias value
    solution  = ChargeTransport.solve(ctsys, inival = solEQ, times=(0.0, tEndQuasiStatic), control = control)
    inival    = solution[end]; initialCond = inival

    ## here we calculate the quasi static solution
    solution2 = ChargeTransport.solve(ctsys, inival = inival, times=(tEndQuasiStatic, 250), control = control)
    solStead  = solution2[end]

    ################################################################################
    if test == false
        println("Calculate error and relative entropy")
    end
    ################################################################################

    tvalues = solution2.t; number_tsteps = length(tvalues)
    p       = data.params

    #######################
    # save relative entropy w.r.t. steady state
    relEntropy  = zeros(0); relEntropyPsi = zeros(0); relEntropyn = zeros(0)
    relEntropyp = zeros(0); relEntropya   = zeros(0)

    # save L2 errors
    totall2Err  = zeros(0); psil2Err      = zeros(0); iphinl2Err  = zeros(0)
    iphipl2Err  = zeros(0); iphial2Err    = zeros(0)

    for istep = 1:number_tsteps

        sol = solution2[istep]

        E_psi   = 0.0; E_n   = 0.0; E_p   = 0.0; E_a   = 0.0 # entropy w.r.t. steady state
        err_psi = 0.0; err_n = 0.0; err_p = 0.0; err_a = 0.0 # L2 error

        #### relative entropy calculations ####
        inode   = 1
        Nn      = data.params.densityOfStates[iphin, regionflake]
        Np      = data.params.densityOfStates[iphip, regionflake]
        Na      = data.params.densityOfStates[iphia, regionflake]

        # calculate first values
        n1      = get_density(sol,      data, iphin, regionflake, inode = inode)
        p1      = get_density(sol,      data, iphip, regionflake, inode = inode)
        a1      = get_density(sol,      data, iphia, regionflake, inode = inode)
        nStead1 = get_density(solStead, data, iphin, regionflake, inode = inode)
        pStead1 = get_density(solStead, data, iphip, regionflake, inode = inode)
        aStead1 = get_density(solStead, data, iphia, regionflake, inode = inode)

        for ix = 1:length(coord)-1

            h            = coord[ix+1] - coord[ix]
            gradpsiStead = (solStead[ipsi, ix+1] - solStead[ipsi, ix]) / h
            gradpsi      = (sol[ipsi, ix+1] - sol[ipsi, ix]) / h

            n2      = get_density(sol,      data, iphin, regionflake, inode = ix)
            p2      = get_density(sol,      data, iphip, regionflake, inode = ix)
            a2      = get_density(sol,      data, iphia, regionflake, inode = ix)
            nStead2 = get_density(solStead, data, iphin, regionflake, inode = ix)
            pStead2 = get_density(solStead, data, iphip, regionflake, inode = ix)
            aStead2 = get_density(solStead, data, iphia, regionflake, inode = ix)

            # evaluation of entropy contributions
            E_psi = E_psi + 0.5 * data.params.dielectricConstant[regionflake] * h * ( gradpsi - gradpsiStead )^2
            E_n   = E_n   + 0.5 * h * ( HElectric(p, n1, nStead1, iphin, regionflake) + HElectric(p, n2, nStead2, iphin, regionflake) )
            E_p   = E_p   + 0.5 * h * ( HElectric(p, p1, pStead1, iphip, regionflake) + HElectric(p, p2, pStead2, iphip, regionflake) )
            E_a   = E_a   + 0.5 * h * ( HIonic(p, a1, aStead1, iphia, regionflake) + HIonic(p, a2, aStead2, iphia, regionflake) )

            # evaluation of L2 error contributions
            err_psi = err_psi + h * abs(sol[ipsi,  ix] - solStead[ipsi,  ix])^2
            err_n   = err_n   + h * abs(sol[iphin, ix] - solStead[iphin, ix])^2
            err_p   = err_p   + h * abs(sol[iphip, ix] - solStead[iphip, ix])^2
            err_a   = err_a   + h * abs(sol[iphia, ix] - solStead[iphia, ix])^2

            n1 = n2; nStead1 = nStead2
            p1 = p2; pStead1 = pStead2
            a1 = a2; aStead1 = aStead2

        end

        push!(relEntropyPsi, E_psi); push!(relEntropyn, E_n); push!(relEntropyp, E_p)
        push!(relEntropya,   E_a);   push!(relEntropy,  E_psi + E_n + E_p + E_a)

        err_tot = err_psi + err_n + err_p + err_a

        push!(totall2Err, err_tot); push!(psil2Err,   err_psi); push!(iphinl2Err, err_n)
        push!(iphipl2Err, err_p);   push!(iphial2Err, err_a)

    end

    naInt = get_density(solution2[end], data, iphia, regionflake, inode = 1:length(coord))
    naAvg = trapz(coord,  naInt)./h_flake

    println("Average vacancy concentration: $naAvg 1/m^3");

    ##############################################################

    function adjustWithTolerance!(xx)
        eps = 5.0e-18
        for i in eachindex(xx)
            if xx[i] <= eps
                xx[i] = eps
            end
        end
        return xx
    end

    relEntropy  = adjustWithTolerance!(relEntropy); psil2Err2   = adjustWithTolerance!(psil2Err)
    iphinl2Err2 = adjustWithTolerance!(iphinl2Err); iphipl2Err2 = adjustWithTolerance!(iphinl2Err);
    iphial2Err2 = adjustWithTolerance!(iphinl2Err)

    if naAvgSweep
        try

            trelEntrop = 0.0; tpsiErr  = 0.0; tphinErr = 0.0
            tphipErr   = 0.0; tphiaErr = 0.0

            trelEntrop = tvalues[findall(x->x==5.0e-18, relEntropy)[1]];  tpsiErr  = tvalues[findall(x->x==5.0e-18, psil2Err2)[1]]
            tphinErr   = tvalues[findall(x->x==5.0e-18, iphinl2Err2)[1]]; tphipErr = tvalues[findall(x->x==5.0e-18, iphipl2Err2)[1]]
            tphiaErr   = tvalues[findall(x->x==5.0e-18, iphial2Err2)[1]]

            println("For rel entropy, convergence time: $trelEntrop s"); println("For psil2Err,   convergence time: $tpsiErr s")
            println("For iphinl2Err,  convergence time: $tphinErr s");   println("For iphipl2Err, convergence time: $tphipErr s")
            println("For iphial2Err,  convergence time: $tphiaErr s")

            return [trelEntrop, tpsiErr, tphinErr, tphipErr, tphiaErr], naAvg

        catch y
                @warn("Exception: ", y) # only warn and continue in case of an error
        end

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

        nn   = get_density(solStead, data, iphin, regionflake, inode = 1:length(coord))
        np   = get_density(solStead, data, iphip, regionflake, inode = 1:length(coord))
        na   = get_density(solStead, data, iphia, regionflake, inode = 1:length(coord))

        nnIC = get_density(initialCond, data, iphin, regionflake, inode = 1:length(coord))
        npIC = get_density(initialCond, data, iphip, regionflake, inode = 1:length(coord))
        naIC = get_density(initialCond, data, iphia, regionflake, inode = 1:length(coord))

        subplot(311)
        PyPlot.semilogy(coord./μm, nn, color= "green", linewidth = 3, label =  "\$ n_{\\mathrm{n}}^\\infty \$" )
        PyPlot.semilogy(coord./μm, nnIC, linewidth = 3, linestyle= ":", color="black", label= "IC" )
        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=16)

        PyPlot.ylim(5.0e20, 5.0e28)
        ##################################
        subplot(312)
        PyPlot.semilogy(coord./μm, np, color= "red",   linewidth = 3, label =  "\$ n_{\\mathrm{p}}^\\infty\$" )
        PyPlot.semilogy(coord./μm, npIC, linewidth = 3, linestyle= ":", color="black")
        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=16)

        PyPlot.ylim(1.0e-12, 1.0e10)
        ##################################
        subplot(313)
        PyPlot.semilogy(coord./μm, na, color= "gold",  linewidth = 3, label =  "\$ n_{\\mathrm{a}}^\\infty\$" )
        PyPlot.semilogy(coord./μm, naIC, linewidth = 3, linestyle= ":", color="black")
        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=16)

        PyPlot.ylim(1.0e-4, 1.0e32)
        PyPlot.yticks( [1.0e-4, 1.0e23, 1.0e32] )
        PyPlot.tight_layout()

        if saveFig
            savefig("steady-state-TMDC-dens-Ea-$(Ea2/q).pdf")
        end

        ##################################
        PyPlot.figure() ## Plot initial condition and steady state potentials

        subplot(311)
        PyPlot.plot(coord./μm, solStead[iphin, :], color= "green",  linewidth = 3, label =  "\$ \\varphi_{\\mathrm{n}}^{\\infty}\$" )
        PyPlot.plot(coord./μm, solStead[iphip, :], color= "red",   linewidth = 3, label =  "\$ \\varphi_{\\mathrm{p}}^{\\infty}\$" )

        PyPlot.plot(coord./μm, initialCond[iphin, :], linewidth = 3, linestyle= ":", color="black", label= "IC" )
        PyPlot.plot(coord./μm, initialCond[iphip, :], linewidth = 3, linestyle= ":", color="black")

        PyPlot.grid()
        PyPlot.ylim(-0.15, 1.75)
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("potential [V]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=14)
        PyPlot.tight_layout()

        ####################

        subplot(312)
        PyPlot.plot(coord./μm, solStead[iphia, :], color= "gold", linewidth = 3, label =  "\$ \\varphi_{\\mathrm{a}}^{\\infty}\$" )
        PyPlot.plot(coord./μm, initialCond[iphia, :], linewidth = 3, linestyle= ":", color="black" )

        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("potential [V]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=14)
        PyPlot.ylim(-2.0, 2.1) #-0.55
        PyPlot.tight_layout()

        ####################

        subplot(313)
        PyPlot.plot(coord./μm, solStead[ipsi,:], label = "\$\\psi^{\\infty}\$", color="b", linewidth= 3)
        PyPlot.plot(coord./μm, initialCond[ipsi, :], linewidth = 3, linestyle= ":", color="black")

        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("potential [V]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=16)
        PyPlot.ylim(-4.3, -1.4)
        #PyPlot.ylim(-4.3, -2.25)
        PyPlot.tight_layout()
        if saveFig
            savefig("steady-state-TMDC-sol-Ea-$(Ea2/q).pdf")
        end

    end

    if plotting

        ################################################################################
        #  Putting a tolerance
        ################################################################################

        PyPlot.figure()

        itime = 3:10:length(tvalues)

        PyPlot.semilogy(tvalues[itime], relEntropy[itime],  color = "darkgreen", marker = "x")
        PyPlot.xlabel("time [s]", fontsize=17)
        PyPlot.ylabel("relative entropy w.r.t. steady state", fontsize=17)
        PyPlot.grid()
        PyPlot.tight_layout()
        PyPlot.ylim(0.2e-18, 1.0e3)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()

        if saveFig
            savefig("rel-entropy-TMDC-Ea-$(Ea2/q).pdf")
        end

        ######################################################

        PyPlot.figure()

        PyPlot.semilogy(tvalues[itime], psil2Err2[itime],   color = "blue", linewidth =4, marker = "o", label = "\$ \\psi  - \\psi^\\infty \$")
        PyPlot.semilogy(tvalues[itime], iphinl2Err2[itime], color = "green",  linewidth =4,  label = "\$ \\varphi_{\\mathrm{n}} - \\varphi_{\\mathrm{n}}^\\infty \$")
        PyPlot.semilogy(tvalues[itime], iphipl2Err2[itime], color = "red", linewidth =4,  linestyle=":",  label = "\$ \\varphi_{\\mathrm{p}} - \\varphi_{\\mathrm{p}}^\\infty \$")
        PyPlot.semilogy(tvalues[itime], iphial2Err2[itime], color = "gold", linestyle="dashed", marker = "o", label = "\$ \\varphi_{\\mathrm{a}} - \\varphi_{\\mathrm{a}}^\\infty \$")
        PyPlot.xlabel("time [s]", fontsize=17)
        PyPlot.ylabel("\$ || u - u^{\\infty} ||_{L^2}^2\$", fontsize=17)
        PyPlot.legend(fontsize=16)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.grid()
        PyPlot.tight_layout()
        PyPlot.ylim(0.2e-18, 1.0e3)
        ##########################

        if saveFig
            savefig("l2-error-TMDC-Ea-$(Ea2/q).pdf")
        end

    end

    testval = sum(filter(!isnan, solStead))/length(solStead) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
      testval = -0.6240571264613662
    main(test = true) ≈ testval
end


end # module
