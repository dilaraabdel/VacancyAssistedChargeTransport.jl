
#=
# Large time behavior of a drift-diffusion model
ressembling a perovskite model with all parameters to one or zero.
Boundary conditions are constant.

=#

module TestCase1A

using ChargeTransport
using ExtendableGrids
using PyPlot
using LinearAlgebra

ChargeTransport.set_unity_constants()

include("../parameters/Params_unity.jl")

function main(;plotting = false, saveFig = false, test = false)

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    ################################################################################
    if test == false
        println("Define physical parameters and model")
    end
    ################################################################################

    # doping
    Cn             =  0.1
    Cp             = -0.1
    Ca             = -0.1
    # contact voltage
    contactVoltage = 0.5

    # For all other parameters, we refer to the parameters template "Params_unity.jl"

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Set up grid and regions")
    end
    ################################################################################

    n         = 2 * (2^(nref-1))
    hp        = h_pdoping/convert(Float64, n)
    hi        = h_intrinsic/convert(Float64, n)
    hn        = h_ndoping/convert(Float64, n)

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
        println("Define System and fill in information about model")
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
    params.doping[iphip, regionAcceptor]           = Cp
    params.doping[iphia, regionIntrinsic]          = Ca
    params.doping[iphin, regionDonor]              = Cn

    data.params = params
    ctsys       = System(grid, data, unknown_storage=:sparse)
    ipsi        = data.index_psi

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Define initial condition")
    end
    ################################################################################

    # initialize solution and starting vectors
    inival  = unknowns(ctsys); initialCond  = unknowns(ctsys)
    sol     = unknowns(ctsys); solStat      = unknowns(ctsys)
    inival .= 0.0;             initialCond .= 0.0

    ## initial condition for time solver
    initialCond[iphin, :]                   = 0.2 .* sin.(1/ coord[end] .* pi .* coord) .+ 0.5
    initialCond[iphip, :]                   = 0.2 .* sin.(1/ coord[end] .* pi .* coord) .+ 0.5
    initialCond[ipsi, :]                    = 0.2 .* sin.(1/ coord[end] .* pi .* coord) .+ (0.5 + asinh(Cn/2))
    initialCond[iphia, length_p:length_pi] .= 0.1

    inival                                  = initialCond

    ## since the constant which represents the constant quasi Fermi potential of anion vacancies is undetermined, we need
    ## to fix it in the following
    ctsys.fvmsys.boundary_factors[iphia, bregionJ2] = 1.0e30
    ctsys.fvmsys.boundary_values[iphia, bregionJ2]  = 0.0

    ## We use the equilibrium solution as starting guess for the stationary solution
    solEQ                                           = solve(ctsys, inival = inival)

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Calculate stationary solution")
    end
    ################################################################################

    ## set calculation type to outOfEquilibrium calculations
    data.calculationType = OutOfEquilibrium

    ## Define boundary conditions
    set_contact!(ctsys, bregionAcceptor, Δu = contactVoltage)
    set_contact!(ctsys, bregionDonor,    Δu = contactVoltage)

    ## We do not make use of the stationary solution in the following, it is just used for a consistency check at the end
    ## to see, if steady state and stationary solution coincide
    solStat = solve(ctsys, inival = solEQ)

    ## save the Dirichlet boundary values
    psiD = solStat[ipsi, 1];                  phiD = solStat[iphin, 1]
    nD   = data.F[iphin](zn * (phiD - psiD)); pD   = data.F[iphip](zp * (phiD - psiD))

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Time Loop")
    end
    ################################################################################

    ## initialize Array saving solutions for each time step
    solphin       = zeros(Float64, length(coord), numbertsteps)
    solphip       = zeros(Float64, length(coord), numbertsteps)
    solphia       = zeros(Float64, length(coord), numbertsteps)
    solpsi        = zeros(Float64, length(coord), numbertsteps)

    solphin[:, 1] = inival[iphin, :]; solphip[:, 1] = inival[iphip, :]
    solphia[:, 1] = inival[iphia, :]; solpsi[ :, 1] = inival[ipsi,  :]

    for istep = 2:numbertsteps

        t  = tvalues[istep]       # Actual time
        Δt = t - tvalues[istep-1] # Time step size

        if test == false
            println("          time value: Δt = $(t)")
        end

        sol    = solve(ctsys; inival = inival, tstep = Δt)
        inival = sol

        solphin[:, istep] .= sol[iphin, :]; solphip[:, istep] .= sol[iphip, :]
        solphia[:, istep] .= sol[iphia, :]; solpsi[ :, istep] .= sol[ipsi,  :]

    end

    if test == false
        println("*** done\n")
    end

    ################################################################################
    if test == false
        println("Calculate error and relative entropy")
    end
    ################################################################################

    regionNodes  = [2:length_p, length_p+1:length_pi, length_pi+1:length(coord)-1] # first node will be calculated seperately

    #######################
    # save relative entropy w.r.t. steady state
    relEntropy   = zeros(0); relEntropyPsi  = zeros(0); relEntropyn  = zeros(0)
    relEntropyp  = zeros(0); relEntropya    = zeros(0)

    # save relative entropy w.r.t. BC
    relEntropyD  = zeros(0); relEntropyPsiD = zeros(0); relEntropynD = zeros(0)
    relEntropypD = zeros(0); relEntropyaD   = zeros(0)

    # save L2 errors
    totall2Err   = zeros(0); psil2Err       = zeros(0); iphinl2Err   = zeros(0)
    iphipl2Err   = zeros(0); iphial2Err     = zeros(0)

    solStead     = [solphin[:, end] solphip[:, end] solphia[:, end] solpsi[:, end]]'

    for istep = 1:length(tvalues)

        sol    = [solphin[:, istep] solphip[:, istep] solphia[:, istep] solpsi[:, istep]]'

        E_psi  = 0.0; E_n  = 0.0; E_p  = 0.0; E_a  = 0.0 # entropy w.r.t. steady state
        E_psiD = 0.0; E_nD = 0.0; E_pD = 0.0; E_aD = 0.0 # entropy w.r.t. BC

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
                gradpsiD     = 0.0
                gradpsi      = (sol[ipsi, ix+1] - sol[ipsi, ix]) / h

                n2           = get_density(sol, data, iphin, ireg, inode = ix)
                p2           = get_density(sol, data, iphip, ireg, inode = ix)
                nStead2      = get_density(solStead, data, iphin, ireg, inode = ix)
                pStead2      = get_density(solStead, data, iphip, ireg, inode = ix)

                # evaluation of entropy contributions w.r.t steady state
                E_psi  = E_psi + 0.5 * h * ( gradpsi - gradpsiStead )^2
                E_n    = E_n   + 0.5 * h * ( H(n1/Nn, nStead1/Nn) + H(n2/Nn, nStead2/Nn) )
                E_p    = E_p   + 0.5 * h * ( H(p1/Np, pStead1/Np) + H(p2/Np, pStead2/Np) )

                # evaluation of entropy contributions w.r.t BC
                E_psiD = E_psiD + 0.5 * h * ( gradpsi - gradpsiD )^2
                E_nD   = E_nD   + 0.5 * h * ( H(n1/Nn, nD/Nn) + H(n2/Nn, nD/Nn) )
                E_pD   = E_pD   + 0.5 * h * ( H(p1/Np, pD/Np) + H(p2/Np, pD/Np) )

                if ireg == regionIntrinsic # only in intrinsic layer
                    a2      = get_density(sol, data, iphia, ireg, inode = ix)
                    aStead2 = get_density(solStead, data, iphia, ireg, inode = ix)

                    E_aD    = E_aD + 0.5 * h * ( primitiveFDM1Inv(a1/Na)   + primitiveFDM1Inv(a2/Na) )
                    E_a     = E_a  + 0.5 * h * ( HFDM1(a1/Na, aStead1/Na) + HFDM1(a2/Na, aStead2/Na) )

                    a1      = a2; aStead1 = aStead2
                else
                    E_aD    = E_aD  + 0.0; E_a = E_a   + 0.0
                end

                n1 = n2; nStead1 = nStead2
                p1 = p2; pStead1 = pStead2

            end

        end

        push!(relEntropyPsi, E_psi); push!(relEntropyn, E_n); push!(relEntropyp, E_p)
        push!(relEntropya,   E_a);   push!(relEntropy,  E_psi + E_n + E_p + E_a)

        push!(relEntropyPsiD, E_psiD); push!(relEntropynD, E_nD); push!(relEntropypD, E_pD)
        push!(relEntropyaD,   E_aD);   push!(relEntropyD,  E_psiD + E_nD + E_pD + E_aD)

        #### L2 error calculation ####
        # we calculate the L2 error on each subdomain
        err_psi =   sqrt(hp) .* norm( sol[ipsi,  1:length_p-1]       - solStead[ipsi,  1:length_p-1],       2)
                  + sqrt(hi) .* norm( sol[ipsi,  length_p:length_pi] - solStead[ipsi,  length_p:length_pi], 2)
                  + sqrt(hn) .* norm( sol[ipsi,  length_pi+1:end]    - solStead[ipsi,  length_pi+1:end],    2)
        ####################
        err_n   =   sqrt(hp) .* norm( sol[iphin, 1:length_p-1]       - solStead[iphin, 1:length_p-1],       2)
                  + sqrt(hi) .* norm( sol[iphin, length_p:length_pi] - solStead[iphin, length_p:length_pi], 2)
                  + sqrt(hn) .* norm( sol[iphin, length_pi+1:end]    - solStead[iphin, length_pi+1:end],    2)
        ####################
        err_p   =   sqrt(hp) .* norm( sol[iphip, 1:length_p-1]       - solStead[iphip, 1:length_p-1],       2)
                  + sqrt(hi) .* norm( sol[iphip, length_p:length_pi] - solStead[iphip, length_p:length_pi], 2)
                  + sqrt(hn) .* norm( sol[iphip, length_pi+1:end]    - solStead[iphip, length_pi+1:end],    2)
        ####################
        err_a   =   sqrt(hi) .* norm( sol[iphia, length_p:length_pi] - solStead[iphia, length_p:length_pi], 2)

        err_tot = err_psi + err_n + err_p + err_a

        push!(totall2Err, err_tot); push!(psil2Err,   err_psi); push!(iphinl2Err, err_n)
        push!(iphipl2Err, err_p);   push!(iphial2Err, err_a)

    end # time loop

    if test == false
        println("*** done\n")
    end

    if plotting
        ################################################################################
        if test == false
            println("Plotting")
        end
        ################################################################################
        PyPlot.figure() ## Plot steady state densities

        # we have for electrons and holes a homogeneous set of parameters
        nn = get_density(sol, data, iphin, regionDonor, inode = 1:length(coord))
        np = get_density(sol, data, iphip, regionDonor, inode = 1:length(coord))
        na = get_density(sol, data, iphia, regionIntrinsic, inode = 1:length(coord))

        PyPlot.plot(coord, nn, color = "green", linewidth = 3, label = "\$ n_{\\mathrm{n}}^{\\infty} \$")
        PyPlot.plot(coord, np, color = "red",   linewidth = 3, label = "\$ n_{\\mathrm{p}}^{\\infty} \$")
        PyPlot.plot(coord, na, color = "gold",  linewidth = 3, label = "\$ n_{\\mathrm{a}}^{\\infty} \$")

        PyPlot.ylim(0.3, 1.3)
        PyPlot.grid()
        PyPlot.xlabel("space", fontsize=17)
        PyPlot.ylabel("density", fontsize=17)
        PyPlot.legend(fontsize=16)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()

        if saveFig
            savefig("steady-state-test-case-1A-dens.pdf")
        end

        ##################################
        PyPlot.figure() ## Plot initial condition and steady state potentials

        subplot(211)
        PyPlot.plot(coord, sol[iphin, :], color= "green", linestyle = "-", linewidth = 3, label =  "\$ \\varphi_{\\mathrm{n}}^{\\infty}\$" )
        PyPlot.plot(coord, sol[iphip, :], color= "red",   linestyle = ":", linewidth = 3, label =  "\$ \\varphi_{\\mathrm{p}}^{\\infty}\$" )
        PyPlot.plot(coord, sol[iphia, :], color= "gold",  linewidth = 3, label =  "\$ \\varphi_{\\mathrm{a}}^{\\infty}\$" )

        PyPlot.plot(coord, initialCond[iphin, :], linewidth = 3, linestyle= ":", color="black", label= "IC" )
        PyPlot.plot(coord[length_p:length_pi], initialCond[iphia, length_p:length_pi], linewidth = 3, linestyle= ":", color="black" )

        PyPlot.grid()
        PyPlot.xlabel("space", fontsize=17)
        PyPlot.ylabel("potential", fontsize=17)
        PyPlot.legend( fontsize=14)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.ylim(-0.1, 0.8)
        PyPlot.tight_layout()

        ########################################
        subplot(212)
        PyPlot.plot(coord, sol[ipsi,:], label = "\$\\psi^{\\infty}\$", color="b", linewidth= 3)
        PyPlot.plot(coord, initialCond[ipsi, :], linewidth = 3, linestyle= ":", color="black", label= "IC")

        PyPlot.ylim(0.5, 0.8)
        PyPlot.grid()
        PyPlot.xlabel("space", fontsize=17)
        PyPlot.ylabel("potential", fontsize=17)
        PyPlot.legend(fontsize=14)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()
        if saveFig
            savefig("steady-state-test-case-1A-sol.pdf")
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

        PyPlot.figure()

        plotStep = 19
        PyPlot.semilogy(tvalues[2:plotStep:end], adjustWithTolerance!(relEntropyD[2:plotStep:end]), color = "darkblue", marker = "o", linewidth =3, linestyle = ":", label = "w.r.t. boundary values")
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!(relEntropy[1:plotStep:end]),  color = "darkgreen", marker = "x", linewidth =3, markersize = 9, linestyle = ":", label = "w.r.t. steady state")
        PyPlot.xlabel("time", fontsize=17)
        PyPlot.ylabel("Relative entropy", fontsize=17)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.grid()
        PyPlot.legend(fontsize=16)
        PyPlot.xlim(-3.0, 83)
        PyPlot.ylim(0.2e-18, 4.0e0)
        PyPlot.tight_layout()
        if saveFig
            savefig("rel-entropy-test-case-1A.pdf")
        end
        ######################################################

        PyPlot.figure()
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!( (psil2Err[1:plotStep:end]).^2 ),  color = "blue",  marker = "x", markersize = 9, linewidth =4, label = "\$ \\psi  - \\psi^\\infty \$")
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!( (iphinl2Err[1:plotStep:end]).^2 ), color = "green", marker = "x", markersize = 9, linewidth =4,  label = "\$ \\varphi_{\\mathrm{n}} - \\varphi_{\\mathrm{n}}^\\infty \$")
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!( (iphipl2Err[1:plotStep:end]).^2 ), color = "red",  linestyle="dashed",  marker = "o", label = "\$ \\varphi_{\\mathrm{p}} - \\varphi_{\\mathrm{p}}^\\infty \$")
        PyPlot.semilogy(tvalues[1:plotStep:end], adjustWithTolerance!( (iphial2Err[1:plotStep:end]).^2 ), color = "gold", linestyle="dashed", marker = "o", label = "\$ \\varphi_{\\mathrm{a}} - \\varphi_{\\mathrm{a}}^\\infty \$")
        PyPlot.xlabel("time", fontsize=17)
        PyPlot.xlim(-3.0, 83)
        PyPlot.ylim(0.2e-18, 4.0e0)
        PyPlot.ylabel("\$ || u - u^{\\infty} ||_{L^2}^2\$", fontsize=17)
        PyPlot.grid()
        PyPlot.legend(fontsize=16)
        PyPlot.tick_params(which ="both", labelsize=18)
        PyPlot.tight_layout()
        if saveFig
            savefig("l2-error-test-case-1A.pdf")
        end

        if test == false
            println("*** done\n")
        end
    end

    testval = sum(filter(!isnan, solStat))/length(solStat) # when using sparse storage, we get NaN values in solution
    return testval

end #  main

function test()
    testval = 0.40103353805618497
    main(test = true) ≈ testval
end


end # module
