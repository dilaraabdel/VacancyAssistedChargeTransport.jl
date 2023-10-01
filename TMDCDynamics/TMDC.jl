

#=

Calculation of I-V characteristics for a TMDC-based memristive device.

=#

module TMDC

using ChargeTransport
using ExtendableGrids
using PyPlot
using LinearAlgebra
using DelimitedFiles
using Roots

function main(;plotting = true,  test = false, parameter_file = "../parameters/Params_TMDC_S1.jl",            # "../parameters/Params_TMDC_S2.jl"
               SteadyStatCalculation  = false, biasValSteadyStat = 1.6 * V,                                   # for steady state calculations
               mobilitySweep          = false, μaSweep           = 5.0e-14 *(m^2)/(V*s),                      # for mobility sweep
               PrintBarriers          = false,                                                                # post process for plotting the barriers
               barrierSweep           = false, barrierLeftSweep = 1.0e-3 *eV, barrierRightSweep = 1.0e-3 *eV, # for barrier sweep
               enableIons             = true)

    include(parameter_file)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    if mobilitySweep == false && barrierSweep == false
        PyPlot.close("all")
    end
    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

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

    # Apply zero voltage on left boundary and a linear scan protocol on right boundary
    contactVoltageFunction = [zeroVoltage, scanProtocol]

    #########################################################################
    if SteadyStatCalculation # in case of specific post processing

        ## we need to adjust the initial guess such that the times are consistent
        ## with the times in the publication
        if biasValSteadyStat == 4.4 * V
            initialGuess = 4.3
        elseif biasValSteadyStat == 9.6 * V
            initialGuess = 4
        elseif biasValSteadyStat < 0.0 * V
            initialGuess = 4
        else
            initialGuess = 0.0
        end

        g(t)            = scanProtocol(t) - biasValSteadyStat
        tEndSteadyStat = find_zero(g, initialGuess) + period

        ## Define scan protocol function
        function scanProtocolQuasiStatic(t)

            if t <= tEndSteadyStat
                biasVal = scanProtocol(t)
            else
                biasVal = biasValSteadyStat
            end

            return biasVal

        end

        contactVoltageFunction = [zeroVoltage, scanProtocolQuasiStatic]

    end

    #########################################################################
    if mobilitySweep
        μa2 = μaSweep
    else
        μa2 = μa
    end

    #########################################################################
    if barrierSweep
        barrierLeft2  = barrierLeftSweep
        barrierRight2 = barrierRightSweep
    else
        barrierLeft2  = barrierLeft
        barrierRight2 = barrierRight
    end

    #########################################################################

    # For all other parameters, we refer to the parameters template

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    # non-uniform grid
    coord1 = geomspace(0.0,       h_flake/2, 5e-4 * h_flake, 2e-2 * h_flake)
    coord2 = geomspace(h_flake/2, h_flake,   2e-2 * h_flake, 5e-4 * h_flake)
    coord  = glue(coord1, coord2)

    grid   = simplexgrid(coord)

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

    if enableIons
        numberOfCarriers2 = numberOfCarriers
    else
        numberOfCarriers2 = 2
    end

    # initialize Data instance and fill in data
    data                            = Data(grid, numberOfCarriers2, contactVoltageFunction = contactVoltageFunction)
    data.modelType                  = Transient
    data.bulkRecombination          = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                             bulk_recomb_Auger = false,
                                                             bulk_recomb_radiative = false,
                                                             bulk_recomb_SRH = false)

    data.boundaryType[bregionLeft]  = boundaryType
    data.boundaryType[bregionRight] = boundaryType
    data.fluxApproximation         .= ExcessChemicalPotential

    if enableIons
        data.F                      = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA, FermiDiracMinusOne]
        enable_ionic_carrier!(data, ionicCarrier = iphia, regions = [regionflake])
    else
        data.F                      = [FermiDiracOneHalfTeSCA, FermiDiracOneHalfTeSCA]
    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define Params and fill in physical parameters")
    end
    ################################################################################

    params                                           = Params(grid, numberOfCarriers2)

    params.temperature                               = T
    params.UT                                        = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                      = zn
    params.chargeNumbers[iphip]                      = zp


    params.dielectricConstant[regionflake]           = εr * ε0
    params.dielectricConstantImageForce[regionflake] = εi * ε0

    ## effective DOS, band-edge energy and mobilities
    params.densityOfStates[iphin, regionflake]       = Nn
    params.densityOfStates[iphip, regionflake]       = Np


    params.bandEdgeEnergy[iphin, regionflake]        = En
    params.bandEdgeEnergy[iphip, regionflake]        = Ep


    params.mobility[iphin, regionflake]              = μn
    params.mobility[iphip, regionflake]              = μp


    if enableIons
        params.chargeNumbers[iphia]                  = za
        params.densityOfStates[iphia, regionflake]   = Na
        params.bandEdgeEnergy[iphia, regionflake]    = Ea
        params.mobility[iphia, regionflake]          = μa2
    end

    # doping
    params.doping[iphin, regionflake]                = Cn


    ## Schottky contact
    params.SchottkyBarrier[bregionLeft]              = barrierLeft2
    params.SchottkyBarrier[bregionRight]             = barrierRight2
    params.bVelocity[iphin, bregionLeft]             = vn
    params.bVelocity[iphin, bregionRight]            = vn
    params.bVelocity[iphip, bregionLeft]             = vp
    params.bVelocity[iphip, bregionRight]            = vp

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
    if test == false
        control.verbose      = "ed"
    else
        control.verbose      = "e"
    end
    control.damp_initial = damp_initial
    control.damp_growth  = damp_growth
    control.max_round    = max_round
    control.maxiters     = maxiters

    control.abstol       = abstol
    control.reltol       = reltol
    control.tol_round    = tol_round

    control.Δu_opt       = Inf
    control.Δt           = Δt
    control.Δt_min       = Δt_min
    control.Δt_max       = Δt_max
    control.Δt_grow      = Δt_grow

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
    sol     = unknowns(ctsys);
    inival .= 0.0;             initialCond .= 0.0

    solEQ   = equilibrium_solve!(ctsys, control = control, nonlinear_steps = 20)
    inival  = solEQ

    if plotting

        nn = get_density(solEQ, data, iphin, regionflake, inode = 1:length(coord))
        np = get_density(solEQ, data, iphip, regionflake, inode = 1:length(coord))

        PyPlot.semilogy(coord./μm, nn, color = "green", linewidth = 3, label = "\$ n_{\\mathrm{n}} \$")
        PyPlot.semilogy(coord./μm, np, color = "red",   linewidth = 3, label = "\$ n_{\\mathrm{p}} \$")
        if enableIons
            na = get_density(solEQ, data, iphia, regionflake, inode = 1:length(coord))
            PyPlot.semilogy(coord./μm, na, color = "gold",  linewidth = 3, label = "\$ n_{\\mathrm{a}} \$")
        end

        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("density [\$\\frac{1}{\\mathrm{m}^3}\$]", fontsize=17)
        PyPlot.legend(fontsize=16)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()
        ######################
        PyPlot.figure()
        PyPlot.plot(coord./μm, solEQ[iphin, :], color= "green", linewidth = 3, label =  "\$ \\varphi_{\\mathrm{n}} \$")
        PyPlot.plot(coord./μm, solEQ[iphip, :], color= "red",   linewidth = 3, label =  "\$ \\varphi_{\\mathrm{p}} \$" )
        PyPlot.plot(coord./μm, solEQ[ipsi, :],  color= "blue",  linewidth = 3, label =  "\$ \\psi \$" )
        if enableIons
            PyPlot.plot(coord./μm, solEQ[iphia, :], color= "gold",  linewidth = 3, linestyle = ":", label =  "\$ \\varphi_{\\mathrm{a}} \$" )
        end

        PyPlot.grid()
        PyPlot.xlabel("space [\$\\mu\$m]", fontsize=17)
        PyPlot.ylabel("potential [V]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=14)
        PyPlot.tight_layout()

    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    if SteadyStatCalculation
        sol = ChargeTransport.solve(ctsys, inival = inival, times=(0.0, tEndSteadyStat), control = control)
    else
        sol = ChargeTransport.solve(ctsys, inival = inival, times=(0.0, tend),            control = control)
    end

    if test == false
        println("*** done\n")
    end

    if SteadyStatCalculation
        sol2 = ChargeTransport.solve(ctsys, inival = sol[end], times=(tEndSteadyStat, 2*tEndSteadyStat), control = control)
        return coord, sol2, data
    end

    ################################################################################
    #########  IV curve calculation
    ################################################################################

    IV            = zeros(0) # for saving I-V data

    tvalues       = sol.t
    number_tsteps = length(tvalues)
    biasValues    = scanProtocol.(tvalues)

    factory       = TestFunctionFactory(ctsys)
    tf            = testfunction(factory, [bregionLeft], [bregionRight])

    push!(IV, 0.0)
    for istep = 2:number_tsteps
        Δt       = tvalues[istep] - tvalues[istep-1] # Time step size
        inival   = sol[istep-1]
        solution = sol[istep]

        I        = integrate(ctsys, tf, solution, inival, Δt)

        current = 0.0
        for ii = 1:numberOfCarriers2+1
            current = current + I[ii]
        end

        push!(IV, current)

    end

    if PrintBarriers
        return sol, IV, biasValues
    end

    if mobilitySweep || barrierSweep
        itime1 = findall(x->x>=period,   tvalues)[1]
        itime2 = findall(x->x>=2*period, tvalues)[1]

        if mobilitySweep
            writedlm("mobility-sweep-$μaSweep.dat", [biasValues[itime1:itime2] abs.(Area .* IV[itime1:itime2])])
        elseif barrierSweep
            writedlm("barrier-sweep-left-$(barrierLeftSweep/q)-right-$(barrierRightSweep/q).dat", [biasValues[itime1:itime2] abs.(Area .* IV[itime1:itime2])])
        end

        return biasValues[itime1:itime2], abs.(Area .* IV[itime1:itime2])
    end

    if plotting
        PyPlot.figure()
        PyPlot.plot(tvalues, biasValues, marker = "x")
        PyPlot.xlabel("time [s]", fontsize=17)
        PyPlot.ylabel("voltage [V]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.grid()
        PyPlot.tight_layout()

        colors     = ["green", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"]
        linestyles = [":", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-"]

        IV_meas = readdlm("data/Li2018-$measurement.dat")
        PyPlot.figure()

        for xx in 1:length(IV_meas[:, 1])
            PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "4", color = "black")
            if xx == length(IV_meas[:, 1])
                PyPlot.plot(IV_meas[xx, 1], abs.(IV_meas[xx, 2]), marker = "o", markersize = "4", color = "black", label ="measurement")
            end
        end

        # find index of other loop
        for iLoop = 1:numberLoops
            itime1 = findall(x->x>=(iLoop-1) * period, tvalues)[1]
            itime2 = findall(x->x>=iLoop     * period, tvalues)[1]

            PyPlot.semilogy(biasValues[itime1:itime2], abs.(Area .* IV[itime1:itime2]), linewidth = 5, color = colors[iLoop], linestyle = linestyles[iLoop], label ="Cycle $iLoop")

            # writedlm("IV-$paramsname-cycle-$iLoop.dat", [biasValues[itime1:itime2] abs.(Area .* IV[itime1:itime2])])
            # writedlm("sol-$paramsname-cycle-$iLoop.dat", [coord sol[itime1]'])
        end

        PyPlot.grid()
        PyPlot.legend()
        PyPlot.xlabel("bias [V]", fontsize=17)
        PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()
    end

    testval = sum(filter(!isnan, sol[end]))/length(sol[end]) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
      testval = -1.1040054036623108
    main(test = true) ≈ testval
end


end # module
