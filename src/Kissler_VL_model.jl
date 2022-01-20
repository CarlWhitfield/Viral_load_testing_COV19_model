"""
    Kissler_VL_model.jl

Written by Carl Whitfield, University of Manchester, 2021

Only to be included by definitions.jl
"""

#===================== VL trajectory params  =====================#

#====> Kissler et al 2020 <====#
const beta_peak = 0.7
const alpha_peak = 3.5*beta_peak                      #Peak time from LOD distribution params
const beta_decay = 0.3
const alpha_decay_symp = 10.5*beta_decay              #Decay time to LOD distribution params
const alpha_decay_asymp = 6.7*beta_decay
const peakVL_mean = 7.533                             #Peak viral load distribution params
const peakVL_sd = 1.164
const VL_LOD_Kissler = log10(250) - 0.93733/3.60971   #Limit of Detection in Kissler study
const V_0 = 0.5255                  #from Ferretti (initial VL) (if using V0_model_opt)
const onset_frac = 0.38             #frac of peak time added before peak (if using PVT_model_opt)

#Distributions used to generate params
const PVLdist = truncated(Normal(peakVL_mean,peakVL_sd),VL_LOD_Kissler,Inf)
const PTdist = Gamma(alpha_peak,1.0/beta_peak)
const DTsympdist = Gamma(alpha_decay_symp,1.0/beta_decay)
const DTasympdist = Gamma(alpha_decay_asymp,1.0/beta_decay)

""" 
    generate_VL_params_Kissler(Asymp::Bool)

Generate parameters for VL model. 

## Arguments
`Asymp` = Whether individual is asymptomatic

## Returns:
`Float` = peak viral load (log10 copies/ml)
`Float` = peak viral load time (days since infection)
`Float` = growth rate (1/days)
`Float` = decay rate (1/days)
"""
function generate_VL_params_Kissler(Asymp::Bool)
    PVL = rand(PVLdist)  #peak viral load
    PT = rand(PTdist)    #time from LOD to peak
    
    r = log(10)*(PVL - VL_LOD_Kissler)/PT
    
    if onset_opt == V0_model_opt 
        OT = log(10)*(VL_LOD_Kissler - V_0)/r    #time from infection to LOD
    elseif onset_opt == PVT_model_opt 
        OT = onset_frac*PT
    end
    
    if Asymp
        DT = rand(DTasympdist)
    else
        DT = rand(DTsympdist)
    end
    
    tp = OT+PT
    d = log(10)*(PVL - VL_LOD_Kissler)/DT

    return PVL, tp, r, d
end