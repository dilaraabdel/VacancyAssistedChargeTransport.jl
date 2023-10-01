
#####################################################################
############################ parameters ############################


########## charge carriers ##########

const iphin            = 1 # electron quasi Fermi potential
const iphip            = 2 # hole quasi Fermi potential
const iphia            = 3
const numberOfCarriers = 3 # electrons, holes and anion vacancies

########## device geometry ##########

# region numbers
const regionAcceptor   = 1
const regionIntrinsic  = 2
const regionDonor      = 3
const regions          = [regionAcceptor, regionIntrinsic, regionDonor]
const numberOfRegions  = length(regions)

# boundary region numbers
const bregionAcceptor  = 1
const bregionDonor     = 2
const bregionJ1        = 3
const bregionJ2        = 4

## length domains
const h_pdoping        = 2.00
const h_intrinsic      = 2.00
const h_ndoping        = 2.00
const h_total          = h_pdoping + h_intrinsic + h_ndoping

const nref             = 9 # refinement level grid

########## time mesh ##########
const tend             = 80.0
const tvalues          = collect(0.0:1.0e-1:tend)
const numbertsteps     = length(tvalues)

########## physical values ##########

# charge numbers
const zn               = -1
const zp               = 1
const za               = 1

const T                = 1.0 # temperature

# band edge energies
const En               = 0.0
const Ep               = 0.0
const Ea               = 0.0

# effective densities of density of states
const Nn               = 1.0
const Np               = 1.0
const Na               = 1.0

# mobilities
const μn               = 1.0
const μp               = 1.0
const μa               = 1.0

# relative dielectric permittivity
const ε                = 1.0

#####################################################################
##################### Relative entropy functions ####################

primitiveBoltzInv = (x)              -> x * log(x) - x + 1
primitiveFDM1Inv  = (x)              -> x * log(x) + (1-x) * log(1-x) + log(2)
FDM1Inv           = (x)              -> log(x) - log(1-x)
H                 = (icc, icc_stead) -> primitiveBoltzInv(icc) - primitiveBoltzInv(icc_stead) - log(icc_stead) * (icc - icc_stead)
HFDM1             = (icc, icc_stead) -> primitiveFDM1Inv(icc)  - primitiveFDM1Inv(icc_stead)  - FDM1Inv(icc_stead) * (icc - icc_stead)
