# Parameters from Driftfusion:
# https://github.com/barnesgroupICL/Driftfusion/blob/master/Input_files/pedotpss_mapi_pcbm.csv.


# Default parameters of Ionmonger representing TiO2 | MAPI | spiro-OMeTAD


const paramsname         = "DD"
#####################################################################
############################ parameters ############################

########## charge carriers ##########

const iphin              = 1 # electron quasi Fermi potential
const iphip              = 2 # hole quasi Fermi potential
const iphia              = 3
const numberOfCarriers   = 3 # electrons, holes and anion vacancies

########## device geometry ##########

# region numbers
const regionDonor        = 1
const regionIntrinsic    = 2
const regionAcceptor     = 3
const regions            = [regionDonor, regionIntrinsic, regionAcceptor]
const numberOfRegions    = length(regions)

# boundary region numbers
const bregionDonor       = 1
const bregionAcceptor    = 2
const bregionJ1          = 3
const bregionJ2          = 4

## length domains
const h_ndoping          = 85.0  * nm
const h_intrinsic        = 300.0 * nm
const h_pdoping          = 30.0  * nm
const h_total            = h_ndoping + h_intrinsic + h_pdoping

const nref               = 7 # refinement level grid
const n                  = 2 * (2^(nref-1))
const hn                 = h_ndoping/convert(Float64, n)
const hi                 = h_intrinsic/convert(Float64, n)
const hp                 = h_pdoping/convert(Float64, n)

########## time mesh ##########
# we need the smaller steps since we get convergence issues due to photogeneration
tend                     = 1.5e3 * s
t1                       = range(0.0, 2.0e-3, length = 21)
t2                       = collect(2.0e-3:5.0e-0:tend)
tvalues                  = vcat(collect(t1)[1:end-1], t2)
numbertsteps             = length(tvalues)

########## physical values ##########

## charge numbers
const zn                 = -1
const zp                 = 1
const za                 = 1

## temperature
const T                  = 300.0                          *  K

## band edge energies
const En                 = [-3.8, -3.8,  -3.0]           .*  eV
const Ep                 = [-6.2, -5.4,  -5.1]           .*  eV
const Ea                 = [0.0, -4.5,  0.0]             .*  eV
Ea_i                     = Ea[regionIntrinsic]

## effective densities of density of states
const Nn                 = [1.0e25, 1.0e25,  1.0e26]     ./ (m^3)
const Np                 = [1.0e25, 1.0e25,  1.0e26]     ./ (m^3)
const Na                 = [0.0,    1.21e28, 0.0]        ./ (m^3)
const Na_i               = Na[regionIntrinsic]

## mobilities
const μn                 = [1.0e-7, 2.0e-3, 1.0e-5]      .* (m^2) / (V * s)
const μp                 = [1.0e-7, 2.0e-3, 1.0e-5]      .* (m^2) / (V * s)
const μa                 = [0.0, 1.0e-14,  0.0]          .* (m^2) / (V * s)

## relative dielectric permittivity
const ε                  = [3.0, 23.0, 4.0]              .* 1.0

## radiative recombination
r0                       = [0.0, 3.6e-18, 0.0]           .* cm^3 / s

## life times and trap densities
τn                       = [1.0e100, 1.0e-7, 1.0e100]    .* s
τp                       = [1.0e100, 1.0e-7, 1.0e100]    .* s

## SRH trap energies
const EI                 = [-5.0, -4.60, -4.05]          .* eV

## generation
const incidentPhotonFlux = [0.0, 1.4e21, 0.0]            ./ (m^2 * s)
const absorption         = [0.0, 4.2e6, 0.0]             ./ m

## doping
const Cn                 = 2.09e24                        / (m^3)
const Cp                 = 2.09e24                        / (m^3)
const Ca                 = 1.0e24                         / (m^3)

## contact voltage
contactVoltage           = 1.16                           * V

## bias range
const bias               = range(0, stop = contactVoltage, length = 21)

## auxiliary values for IC
const phiaOC             = 0.2                            * V
const Δψ                 = 3.81                           * V


#####################################################################
#####################        Varying eps         ####################
function vacancy_energy(eps)

    if eps == 0.01
        return -4.66    *  eV
    elseif eps == 0.05
        return -4.61    *  eV
    elseif eps == 0.1
        return -4.59    *  eV
    elseif eps == 0.2
        return -4.56    *  eV
    elseif eps == 0.3
        return -4.53    *  eV
    elseif eps == 0.4
        return -4.51    *  eV
    elseif eps == 0.5
        return -4.477   *  eV
    elseif eps == 0.6
        return -4.44    *  eV
    elseif eps == 0.7
        return -4.38    *  eV
    elseif eps == 0.75
        return  -4.35   *  eV
    elseif eps == 0.8
        return -4.30    *  eV
    elseif eps == 0.9
        return -4.16    *  eV
    elseif eps == 0.95
        return -4.04 * eV
    elseif eps >= 0.999
        return -2.0 * eV
    else
        error("sorry, you need adjust the vacancy energy level manually!")
    end

end

function check_IMData_availability(eps)

    Set = [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9]

    if eps in Set
        return true
    elseif eps >= 0.99
        return true
    else
        return false
    end

end


#####################################################################
##################### Relative entropy functions ####################

const kBT      = kB * T
PrimElectric   = (p, nicc, icc, ireg) -> kBT * nicc * ( log(nicc/p.densityOfStates[icc, ireg]) - 1 ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg] * nicc
PrimElectPrime = (p, nicc, icc, ireg) -> kBT * log(nicc/p.densityOfStates[icc, ireg]) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg]

PrimIonic      = (p, nicc, icc, ireg) -> kBT * ( nicc * log(nicc/p.densityOfStates[icc, ireg]) + (p.densityOfStates[icc, ireg] - nicc) * log(1 - nicc/p.densityOfStates[icc, ireg]) ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg] * nicc

PrimIonicPrime = (p, nicc, icc, ireg) -> kBT * ( log(nicc/p.densityOfStates[icc, ireg]) -  log(1 - nicc/p.densityOfStates[icc, ireg]) ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg]

HElectric      = (p, nicc, niccStead, icc, ireg) -> PrimElectric(p, nicc, icc, ireg) - PrimElectric(p, niccStead, icc, ireg) - PrimElectPrime(p, niccStead, icc, ireg) * (nicc - niccStead)
HIonic         = (p, nicc, niccStead, icc, ireg) -> PrimIonic(p, nicc, icc, ireg)  - PrimIonic(p, niccStead, icc, ireg)   - PrimIonicPrime(p, niccStead, icc, ireg) * (nicc - niccStead)