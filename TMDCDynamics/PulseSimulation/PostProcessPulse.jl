

#=
Plotting the pulse simulation data
=#

module PostProcessPulse

using ChargeTransport
using PyPlot
using DelimitedFiles

include("PulseSimulation.jl")

function main(;test = false, runSim = false, simulation = "Simulation1")

    PyPlot.rc("font", family="serif", size=14)
    PyPlot.rc("mathtext", fontset="dejavuserif")
    PyPlot.close("all")

    include("Params_Pulse_$simulation.jl")

    if runSim
        tvalues, biasValues, IV = PulseSimulation.main(;plotting = true,  test = false, simulation = simulation, runSim = false)

        IV = abs.(Area .* IV)

    else
        tvalues    = readdlm("data-pulse/$simulation-t.dat")
        biasValues = readdlm("data-pulse/$simulation-V.dat")
        IV         = readdlm("data-pulse/$simulation-I.dat")
    end

    ################################################################################
    if test == false
        println("Some plotting")
    end
    ################################################################################

    ### bias value set pulse
    function biasValueSetPulse(t)

        t0      = t % periodSetPulse
        tRead   = 0.0
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
            tRead   = t
        elseif  t0 < tRest + tSetPulse + 2*tRampSetPulse + tRest + 2*tRampReadPulse + tReadPulse
            biasVal = ΔuReadPulse - ΔuReadPulse/tRampReadPulse * (t0 - tRest - tSetPulse - 2*tRampSetPulse - tRest - tRampReadPulse - tReadPulse) # ramp down
        end

        return biasVal, tRead

    end


    ### bias value reset pulse
    function biasValueResetPulse(t)

        t0      = t % periodResetPulse
        tRead   = 0.0
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
            tRead   = t
        elseif  t0 < tRest + tResetPulse + 2*tRampResetPulse + tRest + 2*tRampReadPulse + tReadPulse
            biasVal = ΔuReadPulse - ΔuReadPulse/tRampReadPulse * (t0 - tRest - tResetPulse - 2*tRampResetPulse - tRest - tRampReadPulse - tReadPulse) # ramp down
        end

        return biasVal, tRead

    end

    ## Define scan protocol function
    function scanProtocol(t)

        t0 = t % (tendResetPulse + tendSetPulse)

        if t0 < tendSetPulse
            biasVal, tRead = biasValueResetPulse(t0)
        else
            biasVal, tRead = biasValueSetPulse(t0)
        end
            return biasVal, tRead

    end

    #####################################################
    #Get the voltage values corresponding to t

    iReadPulse        = findall(x -> x == ΔuReadPulse, biasValues)
    iReadPulseReduced = zeros(Int, 0)

    for ii = 1:length(iReadPulse)-1
        if iReadPulse[ii+1][1] != (iReadPulse[ii][1]+1)
            push!(iReadPulseReduced, iReadPulse[ii][1])
        end

    end

    PyPlot.plot(IV[iReadPulseReduced], linewidth = 5, color = "darkblue")
    PyPlot.grid()
    PyPlot.xlabel("pulse number \$N\$", fontsize=17)
    PyPlot.ylabel("read current \$ I_r \$ [A]", fontsize=17)
    PyPlot.tick_params(which ="both", labelsize=18)
    PyPlot.tight_layout()


end #  main


end # module
