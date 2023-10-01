

#####################################################################
############################ parameters ############################

const measurement      = "FigS1i"
const paramsname       = "S1"


## set indices of the quasi Fermi potentials
const iphin            = 1 # electron quasi Fermi potential
const iphip            = 2 # hole quasi Fermi potential
const iphia            = 3 # anion vacancy quasi Fermi potential

const numberOfCarriers = 3 # electrons, holes and anion vacancies


########## device geometry ##########
# region numbers
const regionflake      = 1

# boundary region numbers
const bregionLeft      = 1
const bregionRight     = 2

## length of the conducting channel
const h_flake          = 1.0                                                 * μm

## Area of electrode
const Area             = 1.5e-13                                             * m^2

boundaryType           = SchottkyContact
########## time mesh ##########

      numberLoops      = 2
const amplitude        = 13.0                                                * V
const scanrate         = 5.0                                                 * V/s
const tend             = numberLoops * 4*amplitude/scanrate
const period           = 4*amplitude/scanrate

########## physical values ##########

## charge numbers
const zn               = -1
const zp               = 1
const za               = 1

## temperature
const T                = 300.0                                               *  K

## band edge energies
const En               = -4.0                                                *  eV
const Ep               = -5.3                                                *  eV
      Ea               = -4.32                                               *  eV

## effective densities of density of states
const Nn               = 2*( 2*pi*0.55*mₑ * kB*T/(Planck_constant^2))^(3/2)  / (m^3)
const Np               = 2*( 2*pi*0.71*mₑ * kB*T/(Planck_constant^2))^(3/2)  / (m^3)
const Na               = 1.0e28                                              / (m^3)

## mobilities
const μn               = 2.5e-4                                              * (m^2) / (V*s)
const μp               = 2.5e-4                                              * (m^2) / (V*s)
      μa               = 5.0e-14                                             * (m^2) / (V*s)

## relative dielectric permittivity
const εr               = 10.0

## image force dielectric permittivity
const εi               = εr

## doping
const Cn               = 1.0e21                                              / (m^3)

## values in case of SchottkyContact
      barrierLeft      = 1.0e-3                                              *  eV
      barrierRight     = 1.0e-3                                              *  eV
const An               = 4 * pi * q * 0.55 * mₑ * kB^2 / Planck_constant^3
const Ap               = 4 * pi * q * 0.71 * mₑ * kB^2 / Planck_constant^3
const vn               = An * T^2 / (q*Nn)
const vp               = Ap * T^2 / (q*Np)

#####################################################################
#####################      Newton Parameter      ####################
const damp_initial     = 0.1
const damp_growth      = 1.61 # >= 1
const max_round        = 5
const maxiters         = 500

const abstol           = 1.0e-8
const reltol           = 1.0e-8
const tol_round        = 1.0e-8

const Δu_opt           = Inf
const Δt               = 1.0e-3
const Δt_min           = 1.0e-6
const Δt_max           = 5.0e-1
const Δt_grow          = 1.01
#####################################################################
##################### Relative entropy functions ####################

const kBT      = kB * T
PrimElectric   = (p, nicc, icc, ireg) -> kBT * nicc * ( log(nicc/p.densityOfStates[icc, ireg]) - 1 ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg] * nicc
PrimElectPrime = (p, nicc, icc, ireg) -> kBT * log(nicc/p.densityOfStates[icc, ireg]) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg]

PrimIonic      = (p, nicc, icc, ireg) -> kBT * ( nicc * log(nicc/p.densityOfStates[icc, ireg]) + (p.densityOfStates[icc, ireg] - nicc) * log(1 - nicc/p.densityOfStates[icc, ireg]) ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg] * nicc

PrimIonicPrime = (p, nicc, icc, ireg) -> kBT * ( log(nicc/p.densityOfStates[icc, ireg]) -  log(1 - nicc/p.densityOfStates[icc, ireg]) ) - p.chargeNumbers[icc] * p.bandEdgeEnergy[icc, ireg]

HElectric      = (p, nicc, niccStead, icc, ireg) -> PrimElectric(p, nicc, icc, ireg) - PrimElectric(p, niccStead, icc, ireg) - PrimElectPrime(p, niccStead, icc, ireg) * (nicc - niccStead)
HIonic         = (p, nicc, niccStead, icc, ireg) -> PrimIonic(p, nicc, icc, ireg)  - PrimIonic(p, niccStead, icc, ireg)   - PrimIonicPrime(p, niccStead, icc, ireg) * (nicc - niccStead)