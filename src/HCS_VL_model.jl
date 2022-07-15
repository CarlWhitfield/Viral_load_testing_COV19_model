"""
    HCS_VL_model.jl

Written by Carl Whitfield, University of Manchester, 2022

Only to be included by definitions.jl
"""

#===================== VL trajectory params  =====================#

const HCS_mu_vals = [0.0219, 1.67, 0.488, -0.593]
const HCS_S2_mat = [0.1176 -0.005215 0.006055 0.003001;
                   -0.005215 0.1578 -0.1498 0.1088;
                    0.006055 -0.1498 0.1504 -0.09819;
                    0.003001 0.1088 -0.09819 0.1023]


#Distributions used to generate params
const HCS_param_dist = MvNormal(HCS_mu_vals, HCS_S2_mat)

""" 
    generate_VL_params_HCS()

Generate parameters for VL model. 

## Arguments
`Asymp` = Whether individual is asymptomatic

## Returns:
`Float` = peak viral load (log10 copies/ml)
`Float` = peak viral load time (days since infection)
`Float` = growth rate (1/days)
`Float` = decay rate (1/days)
"""
function generate_VL_params_HCS()
    params = rand(HCS_param_dist)
    VL0 = params[1]   #in log10 copies/ml
    tp, r, d = exp.(params[2:4])    #log of piecewise linear params
    PVL = VL0 + r*tp   #in log10 copies/ml
    
    return PVL, tp, log(10)*r, log(10)*d   #transform to PE model params
end