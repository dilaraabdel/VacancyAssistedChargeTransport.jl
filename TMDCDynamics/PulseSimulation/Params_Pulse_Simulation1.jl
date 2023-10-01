
include("../../parameters/Params_TMDC_S1.jl")


########## Pulse parameters ##########
const ResetPulseNumber = 500
const SetPulseNumber   = 500
const numberLoopsRep   = 3 # number of loops repitition

const ΔuResetPulse     = -5     * V
const ΔuSetPulse       = 2      * V
const ΔuReadPulse      = 0.2    * V

const tResetPulse      = 2      * ms
const tSetPulse        = 2      * ms
const tReadPulse       = 100.0  * μs
const tRest            = 2      * ms

const tRampResetPulse  = 1.5    * ms
const tRampSetPulse    = 1.5    * ms
const tRampReadPulse   = 20.0   * µs

## total time span of all set pulses together
const tendSetPulse     = SetPulseNumber * (tRest + tSetPulse + 2*tRampSetPulse + tRest + 2*tRampReadPulse + tReadPulse)
## total time span of all reset pulses together
const tendResetPulse   = ResetPulseNumber * (tRest + tResetPulse + 2*tRampResetPulse + tRest + 2*tRampReadPulse + tReadPulse)

const endTime          = (tendSetPulse + tendResetPulse)*numberLoopsRep
const periodSetPulse   = tendSetPulse/SetPulseNumber
const periodResetPulse = tendResetPulse/ResetPulseNumber