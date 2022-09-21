# import Pkg
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("Distributions")
# Pkg.add("Random")
# Pkg.add("StatsBase")
# Pkg.add("SpecialFunctions")
# Pkg.add("Plots")
# Pkg.add("ReadOnlyArrays")

using DataFrames
using CSV
using Distributions
using Random
using StatsBase
using SpecialFunctions
using Plots
using HTTP
using Base
using ReadOnlyArrays

#testing protocol definitions
const PCR_mass_protocol, LFD_mass_protocol, LFD_pattern = 1,2,3
#LFD sens definitions
const SC_data_2021, SC_data_Oct, porton_down, porton_down_p3b = 1,2,3,4
#Infectivity model definitions
const ke_inf_model_no, flat_inf_model_no, linear_inf_model_no = 1,2,3
#VL model definitions
const ke_model_no, kissler_model_no, HCS_model_no = 1,2,3   
#options definitons for kissler model
const V0_model_opt, PVT_model_opt = 1,2   
#options for infectivity model
const marks_peakinf_opt, unif_peakinf_opt = 1,2

#default options
VL_model = ke_model_no         #choice of viral load model (Ke et al is default)
LFD_model = SC_data_2021         #choice of LFD sensitivity model (social care binned 2021 data is default)
PCR_sens_max = 0.83           #max PCR sensitivity (Ferretti et al 2021)
Inf_model = ke_inf_model_no    #infectivity model
onset_opt = V0_model_opt        #IF using kissler VL model, option for how to generate onset time
peak_inf_opt = marks_peakinf_opt   #IF using flat or linear inf model, what function is used for peak 
PCR_TaT_scale = 1.0    #scale PCR delays by this factor
VL_inf_corr = false
const DefaultPasymp = 0.5
Pasymp = DefaultPasymp

#===================== Symptoms onset timing =====================#

const Default_symp_beta = 4.84 / 2.6^2
const Default_symp_alpha = 4.84 * Default_symp_beta
const DefaultSympDayMean = mean(Gamma(Default_symp_alpha,1.0/Default_symp_beta))
symp_beta = Default_symp_beta
symp_alpha = Default_symp_alpha

#=================================================================#

include("Kissler_VL_model.jl")
include("Ke_VL_inf_models.jl")
include("testing_models.jl")
include("Alt_inf_models.jl")
include("HCS_VL_model.jl")

"""
    generate_asymptomatic()

Generate and return a boolean value determining whether
somebody will develop symptoms.

## Returns: 
`Bool` = whether or not the person in asymptomatic (true = 
          asymptomatic, false = symptoms)
"""
function generate_asymptomatic(p_asymp::Float64 = Pasymp)
    #randomly generate if people are asymptomatic
    return (rand() < p_asymp)
end