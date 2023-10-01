

#=

Calculation of pulse characteristics for a TMDC-based memristive device.

=#

module PulseSimulation

using ChargeTransport
using ExtendableGrids
using PyPlot
using LinearAlgebra
using DelimitedFiles
using Roots

function main(;plotting = true,  test = false, simulation = "Simulation1", runSim = false)

    include("Params_Pulse_$simulation.jl")

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")
    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    ### bias value set pulse
    function biasValueSetPulse(t)

        t0      = t % periodSetPulse
        biasVal = 0.0

        ### Rest Time ###
        if      t0 < tRest
            biasVal = 0.0
        ######################################
        ### Set Pulse Part ###
        elseif t0 < tRest + tRampSetPulse
            biasVal = ΔuSetPulse/tRampSetPulse * (t0 - tRest)
        elseif  t0 < tRest + tSetPulse + tRampSetPulse
            biasVal = ΔuSetPulse
        elseif  t0 < tRest + tSetPulse + 2*tRampSetPulse
            biasVal = ΔuSetPulse - ΔuSetPulse/tRampSetPulse * (t0 - tRest - tSetPulse - tRampSetPulse)
        ######################################
        ### Rest Time ###
        elseif  t0 < tRest + tSetPulse + 2*tRampSetPulse + tRest
            biasVal = 0.0
        ######################################
        ### Read Pulse Part ###
        elseif  t0 < tRest + tSetPulse + 2*tRampSetPulse + tRest + tRampReadPulse
            biasVal = ΔuReadPulse/tRampReadPulse * (t0 - tRest - tSetPulse - 2*tRampSetPulse - tRest) # ramp up
        elseif  t0 < tRest + tSetPulse + 2*tRampSetPulse + tRest + tRampReadPulse + tReadPulse
            biasVal = ΔuReadPulse
        elseif  t0 < tRest + tSetPulse + 2*tRampSetPulse + tRest + 2*tRampReadPulse + tReadPulse
            biasVal = ΔuReadPulse - ΔuReadPulse/tRampReadPulse * (t0 - tRest - tSetPulse - 2*tRampSetPulse - tRest - tRampReadPulse - tReadPulse) # ramp down
        end

        return biasVal

    end


    ### bias value reset pulse
    function biasValueResetPulse(t)

        t0      = t % periodResetPulse
        biasVal = 0.0

        ### Rest Time ###
        if      t0 < tRest
            biasVal = 0.0
        ######################################
        ### Set Pulse Part ###
        elseif t0 < tRest + tRampResetPulse
            biasVal = ΔuResetPulse/tRampResetPulse * (t0 - tRest)
        elseif  t0 < tRest + tResetPulse + tRampResetPulse
            biasVal = ΔuResetPulse
        elseif  t0 < tRest + tResetPulse + 2*tRampResetPulse
            biasVal = ΔuResetPulse - ΔuResetPulse/tRampResetPulse * (t0 - tRest - tResetPulse - tRampResetPulse)
        ######################################
        ### Rest Time ###
        elseif  t0 < tRest + tResetPulse + 2*tRampResetPulse + tRest
            biasVal = 0.0
        ######################################
        ### Read Pulse Part ###
        elseif  t0 < tRest + tResetPulse + 2*tRampResetPulse + tRest + tRampReadPulse
            biasVal = ΔuReadPulse/tRampReadPulse * (t0 - tRest - tResetPulse - 2*tRampResetPulse - tRest) # ramp up
        elseif  t0 < tRest + tResetPulse + 2*tRampResetPulse + tRest + tRampReadPulse + tReadPulse
            biasVal = ΔuReadPulse
        elseif  t0 < tRest + tResetPulse + 2*tRampResetPulse + tRest + 2*tRampReadPulse + tReadPulse
            biasVal = ΔuReadPulse - ΔuReadPulse/tRampReadPulse * (t0 - tRest - tResetPulse - 2*tRampResetPulse - tRest - tRampReadPulse - tReadPulse) # ramp down
        end

        return biasVal

    end

    ## Define scan protocol function
    function scanProtocol(t)

        t0 = t % (tendResetPulse + tendSetPulse)

        if t0 < tendResetPulse
            biasVal = biasValueResetPulse(t0)
        else
            biasVal = biasValueSetPulse(t0)
        end
            return biasVal

    end

    contactVoltageFunction = [zeroVoltage, scanProtocol]

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
    coord1       = geomspace(0.0,       h_flake/2, 5e-4 * h_flake, 2e-2 * h_flake)
    coord2       = geomspace(h_flake/2, h_flake,   2e-2 * h_flake, 5e-4 * h_flake)
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

    data.boundaryType[bregionLeft]  = boundaryType
    data.boundaryType[bregionRight] = boundaryType
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

    params                                           = Params(grid, numberOfCarriers)

    params.temperature                               = T
    params.UT                                        = (kB * params.temperature) / q
    params.chargeNumbers[iphin]                      = zn
    params.chargeNumbers[iphip]                      = zp
    params.chargeNumbers[iphia]                      = za

    params.dielectricConstant[regionflake]           = εr * ε0
    params.dielectricConstantImageForce[regionflake] = εi * ε0

    ## effective DOS, band-edge energy and mobilities
    params.densityOfStates[iphin, regionflake]       = Nn
    params.densityOfStates[iphip, regionflake]       = Np
    params.densityOfStates[iphia, regionflake]       = Na

    params.bandEdgeEnergy[iphin, regionflake]        = En
    params.bandEdgeEnergy[iphip, regionflake]        = Ep
    params.bandEdgeEnergy[iphia, regionflake]        = Ea

    params.mobility[iphin, regionflake]              = μn
    params.mobility[iphip, regionflake]              = μp
    params.mobility[iphia, regionflake]              = μa

    # doping
    params.doping[iphin, regionflake]                = Cn


    ## Schottky contact
    params.SchottkyBarrier[bregionLeft]              = barrierLeft
    params.SchottkyBarrier[bregionRight]             = barrierRight
    params.bVelocity[iphin, bregionLeft]             = vn
    params.bVelocity[iphin, bregionRight]            = vn
    params.bVelocity[iphip, bregionLeft]             = vp
    params.bVelocity[iphip, bregionRight]            = vp

    data.params = params
    ctsys       = System(grid, data, unknown_storage=:sparse)

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
    control.Δt           = 1.0e-5
    control.Δt_min       = Δt_min
    control.Δt_max       = 1.0e-4  # so long is the largest pulse
    control.Δt_grow      = 1.008

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

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("IV Measurement loop")
    end
    ################################################################################

    sol = ChargeTransport.solve(ctsys, inival = inival, times=(0.0, endTime),            control = control)

    if test == false
        println("IV Measurement loop")
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
        for ii = 1:numberOfCarriers+1
            current = current + I[ii]
        end

        push!(IV, current)

    end

    writedlm("Simulation1-I.dat", IV)
    writedlm("Simulation1-V.dat", biasValues)
    writedlm("Simulation1-t.dat", tvalues)

    if runSim
        return tvalues, biasValues, IV
    end

    if plotting
        PyPlot.figure()

        PyPlot.plot(tvalues, biasValues, marker = "x")
        PyPlot.grid()
        PyPlot.xlabel("time [s]", fontsize=17)
        PyPlot.ylabel("voltage [V]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()

        ###########################
        PyPlot.figure()

        PyPlot.semilogy(tvalues, abs.(Area .* IV), linewidth = 5)
        PyPlot.grid()
        PyPlot.xlabel("bias [V]", fontsize=17)
        PyPlot.ylabel("current \$ | I | \$ [A]", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()
    end

    testval = sum(filter(!isnan, sol[end]))/length(sol[end]) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

end # module
