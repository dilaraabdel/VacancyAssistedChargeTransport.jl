# Default parameters of Ionmonger representing TiO2 | MAPI | spiro-OMeTAD

const paramsname         = "IM"
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
const h_ndoping          = 100.0 * nm
const h_intrinsic        = 400.0 * nm
const h_pdoping          = 200.0 * nm
const h_total            = h_ndoping + h_intrinsic + h_pdoping

const nref               = 7 # refinement level grid
const n                  = 2 * (2^(nref-1))
const hn                 = h_ndoping/convert(Float64, n)
const hi                 = h_intrinsic/convert(Float64, n)
const hp                 = h_pdoping/convert(Float64, n)

########## time mesh ##########
const tend               = 2.2e2 * s
const tvalues            = collect(0.0:5.0e-1:tend)
const numbertsteps       = length(tvalues)

########## physical values ##########

## charge numbers
const zn                 = -1
const zp                 = 1
const za                 = 1

## temperature
const T                  = 298.0                          *  K

## band edge energies
const En                 = [-4.0, -3.7, -3.4]            .*  eV
const Ep                 = [-5.8, -5.4, -5.1]            .*  eV
const Ea                 = [0.0, -4.45,  0.0]            .*  eV
const Ea_i               = Ea[regionIntrinsic]

## effective densities of density of states
const Nn                 = [5.0e25, 8.1e24, 5.0e25]      ./ (m^3)
const Np                 = [5.0e25, 5.8e24, 5.0e25]      ./ (m^3)
const Na                 = [0.0,    1.0e27, 0.0]         ./ (m^3)
const Na_i               = Na[regionIntrinsic]

## mobilities
const μn                 = [3.89e-4, 6.62e-3, 3.89e-5]   .* (m^2) / (V * s)
const μp                 = [3.89e-4, 6.62e-3, 3.89e-5]   .* (m^2) / (V * s)
const μa                 = [0.0, 3.93e-16, 0.0]          .* (m^2) / (V * s)

## relative dielectric permittivity
const ε                  = [10.0, 24.1, 3.0]             .* 1.0

## radiative recombination
const r0                 = [6.8e-17, 3.6e-18, 6.3e-17]   .* m^3 / s

## life times and trap densities
const τn                 = [1.0e100, 3.0e-9, 1.0e100]    .* s
const τp                 = [1.0e100, 3.0e-7, 1.0e100]    .* s

## SRH trap energies
const EI                 = [-5.0, -4.55, -4.1]           .* eV

## generation
const incidentPhotonFlux = [0.0, 1.4e21, 0.0]            ./ (m^2 * s)
const absorption         = [0.0, 1.3e7, 0.0]             ./ m

## doping
const Cn                 = 1.00e24                        / (m^3)
const Cp                 = 1.00e24                        / (m^3)
const Ca                 = 1.6e25                         / (m^3)

## contact voltage
const contactVoltage     = 0.9                            * V

## bias range
const bias               = range(0, stop = contactVoltage, length = 21)

## auxiliary values for IC
const phiaOC             = 0.05                           * V
const Δψ                 = 4.1                            * V
#####################################################################
##################### Relative entropy functions ####################

const kBT      = kB * T
PrimElectric   = (p, nicc, icc, ireg) -> kBT * nicc * ( log(nicc/p.densityOfStates[icc, ireg]) - 1 ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg] * nicc
PrimElectPrime = (p, nicc, icc, ireg) -> kBT * log(nicc/p.densityOfStates[icc, ireg]) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg]

PrimIonic      = (p, nicc, icc, ireg) -> kBT * ( nicc * log(nicc/p.densityOfStates[icc, ireg]) + (p.densityOfStates[icc, ireg] - nicc) * log(1 - nicc/p.densityOfStates[icc, ireg]) ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg] * nicc

PrimIonicPrime = (p, nicc, icc, ireg) -> kBT * ( log(nicc/p.densityOfStates[icc, ireg]) -  log(1 - nicc/p.densityOfStates[icc, ireg]) ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg]

HElectric      = (p, nicc, niccStead, icc, ireg) -> PrimElectric(p, nicc, icc, ireg) - PrimElectric(p, niccStead, icc, ireg) - PrimElectPrime(p, niccStead, icc, ireg) * (nicc - niccStead)
HIonic         = (p, nicc, niccStead, icc, ireg) -> PrimIonic(p, nicc, icc, ireg)  - PrimIonic(p, niccStead, icc, ireg)   - PrimIonicPrime(p, niccStead, icc, ireg) * (nicc - niccStead)