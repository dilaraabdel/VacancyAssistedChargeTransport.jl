
#=
Convergence study of a drift-diffusion model ressembling a perovskite model with
all parameters to one or zero.
Boundary conditions are non-constant.
Investigate the error at a specific time between reference solution and solutions on finer becoming meshes.
=#

module TestCase1BConvergence

using ChargeTransport
using ExtendableGrids
using PyPlot
using LinearAlgebra

ChargeTransport.set_unity_constants()

include("../parameters/Params_unity.jl")

# compute error with respect to the coarser mesh
function compute_error(sol, refsol, REFrefinement, nref)

    phin = sol[1, :]; refphin = refsol[1, 1:2^(REFrefinement-nref):end]
    phip = sol[2, :]; refphip = refsol[2, 1:2^(REFrefinement-nref):end]
    phia = sol[3, :]; refphia = refsol[3, 1:2^(REFrefinement-nref):end]
    psi  = sol[4, :]; refpsi  = refsol[4, 1:2^(REFrefinement-nref):end]

    errphin = norm(phin - refphin, 2); errphip = norm(phip - refphip, 2)
    errphia = norm(phia - refphia, 2); errpsi  = norm(psi  - refpsi,  2)
    errtot  = norm([refpsi, refphin, refphip, refphia] - [psi, phin, phip, phia], 2)

    return errphin, errphip, errphia, errpsi, errtot

end

function main(;plotting = false, saveFig = false, test = true)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # doping
    Cn             = 0.5
    Cp             = 0.5
    Ca             = 0.5
    # contact voltages
    contactVoltage = 1.0

    # For all other parameters, we refer to the parameters template "Params_unity.jl"

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    REFrefinement = 9
    refsol        = Array{Float64, numberOfCarriers + 1}

    # for saving error information
    HVec          = zeros(0)
    Errpsi        = zeros(0); Errphin  = zeros(0); Errphip = zeros(0);
    Errphia       = zeros(0); Errtotal = zeros(0)

    for nref = REFrefinement:-1:2

        if test == false
            println("nref = ", nref)
        end

        n  = 2 * (2^(nref-1))
        hp = h_pdoping/convert(Float64, n)
        hi = h_intrinsic/convert(Float64, n)
        hn = h_ndoping/convert(Float64, n)

        # the grid spacing is in every subdomain the same
        if nref != REFrefinement
            push!(HVec, hp)
        end

        coord_p   = collect(0.0:hp:h_pdoping)
        length_p  = length(coord_p)
        coord_i   = collect(h_pdoping:hi:h_pdoping+h_intrinsic)
        length_pi = length_p + length(coord_i) - 1
        coord_n   = collect(h_pdoping+h_intrinsic:hn:h_pdoping+h_intrinsic+h_ndoping)
        coord     = glue(coord_p, coord_i); coord = glue(coord, coord_n)
        grid      = simplexgrid(coord)

        # cellmask! for defining the subregions and assigning region number
        cellmask!(grid, [0.0],                      [h_pdoping],               regionAcceptor)
        cellmask!(grid, [h_pdoping],                [h_pdoping+h_intrinsic],   regionIntrinsic)
        cellmask!(grid, [h_pdoping+h_intrinsic],    [h_total],                 regionDonor)

        bfacemask!(grid, [h_pdoping],               [h_pdoping],               bregionJ1)  # first  inner interface
        bfacemask!(grid, [h_pdoping + h_intrinsic], [h_pdoping + h_intrinsic], bregionJ2)  # second inner interface

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Define ChargeTransportSystem and fill in information about model")
        end
        ################################################################################

        # initialize Data instance and fill in data
        data                               = Data(grid, numberOfCarriers)
        data.modelType                     = Transient
        data.F                             = [Boltzmann, Boltzmann, FermiDiracMinusOne]
        data.bulkRecombination             = set_bulk_recombination(;iphin = iphin, iphip = iphip,
                                                                    bulk_recomb_Auger = false,
                                                                    bulk_recomb_radiative = false,
                                                                    bulk_recomb_SRH = false)
        data.boundaryType[bregionAcceptor] = OhmicContact
        data.boundaryType[bregionDonor]    = OhmicContact
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

        params                                  = Params(grid, numberOfCarriers)

        params.temperature                      = T
        params.UT                               = (kB * params.temperature) / q
        params.chargeNumbers[iphin]             = zn
        params.chargeNumbers[iphip]             = zp
        params.chargeNumbers[iphia]             = za

        for ireg in 1:numberOfRegions # interior region data

            params.dielectricConstant[ireg]     = ε * ε0

            # effective dos, band edge energy and mobilities
            params.densityOfStates[iphin, ireg] = Nn
            params.densityOfStates[iphip, ireg] = Np

            params.bandEdgeEnergy[iphin, ireg]  = En
            params.bandEdgeEnergy[iphip, ireg]  = Ep

            params.mobility[iphin, ireg]        = μn
            params.mobility[iphip, ireg]        = μp

        end

        # vacancy parameter
        params.densityOfStates[iphia, regionIntrinsic] = Na
        params.bandEdgeEnergy[iphia,  regionIntrinsic] = Ea
        params.mobility[iphia,        regionIntrinsic] = μa

        # doping
        params.doping[iphin, regionDonor]              = Cn
        params.doping[iphia, regionIntrinsic]          = Ca
        params.doping[iphip, regionAcceptor]           = Cp

        data.params = params
        ctsys       = System(grid, data, unknown_storage=:dense)

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Define initial condition and calculate stationary solution")
        end
        ################################################################################

        data.calculationType = OutOfEquilibrium

        # initialize solution and starting vectors
        inival  = unknowns(ctsys); initialCond  = unknowns(ctsys)
        sol     = unknowns(ctsys); solStat      = unknowns(ctsys)
        inival .= 0.0;             initialCond .= 0.0

        ## initial condition for time solver
        initialCond[iphin, :]                   = - 5/90 .* coord.^2 .+ 15/90 .* coord .+ 1.0
        initialCond[iphip, :]                   = - 5/90 .* coord.^2 .+ 15/90 .* coord .+ 1.0
        initialCond[data.index_psi, :]          = - 0.0333333333 * (coord .- 1.7373323077).^2 .+ 0.8531443234
        initialCond[iphia, length_p:length_pi] .= 0.05

        inival                                  = initialCond

        ## Define boundary conditions
        set_contact!(ctsys, bregionAcceptor, Δu = contactVoltage)
        set_contact!(ctsys, bregionDonor,    Δu = 0.0)
        ## since the constant which represents the constant quasi Fermi potential of anion vacancies is undetermined, we need
        ## to fix it in the following
        ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 1.0e30
        ctsys.fvmsys.boundary_values[iphia, bregionJ2]  = 0.0

        # save steady state solution
        solStat                                         = solve(ctsys, inival = inival)

        if test == false
            println("*** done\n")
        end

        ################################################################################
        if test == false
            println("Time Loop")
        end
        ################################################################################

        for istep = 2:numbertsteps

            t  = tvalues[istep] # Actual time
            Δt = t - tvalues[istep-1] # Time step size

            sol    = solve(ctsys; inival = inival, tstep = Δt)
            inival = sol

        end # time loop

        if nref != REFrefinement
            errphin, errphip, errphia, errpsi, errtot = compute_error(sol, refsol, REFrefinement, nref)
            push!(Errpsi,   errpsi);  push!(Errphin,  errphin); push!(Errphip, errphip)
            push!(Errphia,  errphia); push!(Errtotal, errtot)

        else
            refsol = sol
        end

        if test == false
            println("*** done\n")
        end

    end # loop

    if plotting
        ################################################################################
        if test == false
            println("Plotting")
        end
        ################################################################################

        PyPlot.loglog(HVec, sqrt.(HVec).*Errtotal, linewidth = 3, markersize = 10, label = " \$|| u_{\\mathrm{ref}}  - u||_{L^2} \$", marker = "o")
        PyPlot.loglog(HVec, 1.3e-2*HVec.*HVec  , "k--", linewidth = 3, label = "\$\\mathcal{O}(h^2)\$")

        PyPlot.grid()
        PyPlot.xlabel(" \$ h \$ ", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.legend(fontsize=16)
        PyPlot.tight_layout()
        if saveFig
            savefig("convergence-test-case-1B.pdf")
        end

    end

    testval = sum(filter(!isnan, HVec))/length(HVec) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = 0.14174107142857142
    main(test = true) ≈ testval
end


end # module
