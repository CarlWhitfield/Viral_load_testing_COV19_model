"""
    Kissler_VL_model.jl

Written by Carl Whitfield, University of Manchester, 2021

Only to be included by definitions.jl
"""
#===================== VL trajectory params  =====================#

# #====> Kissler et al 2020 <====#
# const beta_peak = 0.7
# const alpha_peak = 3.5*beta_peak                      #Peak time from LOD distribution params
# const beta_decay = 0.3
# const alpha_decay_symp = 10.5*beta_decay              #Decay time to LOD distribution params
# const alpha_decay_asymp = 6.7*beta_decay
# const peakVL_mean = 7.533                             #Peak viral load distribution params
# const peakVL_sd = 1.164
# const VL_LOD_Kissler = log10(250) - 0.93733/3.60971   #Limit of Detection in Kissler study
# const V_0 = 0.5255                  #from Ferretti (initial VL) (if using V0_model_opt)
# const onset_frac = 0.38             #frac of peak time added before peak (if using PVT_model_opt)

# #Distributions used to generate params
# const PVLdist = truncated(Normal(peakVL_mean,peakVL_sd),VL_LOD_Kissler,Inf)
# const PTdist = Gamma(alpha_peak,1.0/beta_peak)
# const DTsympdist = Gamma(alpha_decay_symp,1.0/beta_decay)
# const DTasympdist = Gamma(alpha_decay_asymp,1.0/beta_decay)

let filepath = @__FILE__
    current_dir = dirname(filepath)
    Kissler_symp_path = joinpath(current_dir,"df_symp_kissler.csv")
    Kissler_asymp_path = joinpath(current_dir,"df_asymp_kissler.csv")
    if isfile(Kissler_symp_path) && isfile(Kissler_asymp_path)
        print("Loading Kissler dataset...\n")
        global df_symp_kissler = DataFrame(CSV.File(Kissler_symp_path))
        global df_asymp_kissler = DataFrame(CSV.File(Kissler_asymp_path))
    else
        print("Downloading Kissler dataset...\n")
        data_url = "https://raw.githubusercontent.com/gradlab/CtTrajectories/main/output/params_df_split.csv"
        res = HTTP.get(data_url);
        kissler_df = DataFrame(CSV.File(res.body))  #import data
        #covert to PE parameters used here
        VL_lod = (40.93733 - 40)/3.60971 + log10(250);
        VLs = VL_lod .+ kissler_df[:,"dp"]./3.60971
        r_slopes = (VLs .- VL_lod)./kissler_df[:,"wp"];
        d_slopes = (VLs .- VL_lod)./kissler_df[:,"wr"];

        #output symp and asymp data separately
        symp_dict = Dict()
        bool_symp = (kissler_df[:,"symptomatic"].==1)
        symp_dict[:VL] = VLs[bool_symp]
        symp_dict[:gr] = r_slopes[bool_symp]
        symp_dict[:dr] = d_slopes[bool_symp];

        asymp_dict = Dict()
        bool_asymp = (kissler_df[:,"symptomatic"].==0)
        asymp_dict[:VL] = VLs[bool_asymp]
        asymp_dict[:gr] = r_slopes[bool_asymp]
        asymp_dict[:dr] = d_slopes[bool_asymp];

        global df_symp_kissler = DataFrame(symp_dict)
        global df_asymp_kissler = DataFrame(asymp_dict)
        
        CSV.write(Kissler_symp_path, df_symp_kissler)
        CSV.write(Kissler_asymp_path, df_asymp_kissler)
    end
end

const V_0 = 0.5255                  #from Ferretti (initial VL) (if using V0_model_opt)
const onset_frac = 0.38             #frac of peak time added before peak (if using PVT_model_opt)
const VL_LOD_Kissler = log10(250) - 0.93733/3.60971  #Limit of Detection in Kissler study

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
    if Asymp
        df = df_asymp_kissler
    else
        df = df_symp_kissler
    end
    nh = rand(1:size(df,1))
    PVL = df[nh,:VL]
    gr = df[nh,:gr]
    dr = df[nh,:dr]
    if onset_opt == V0_model_opt 
        tp = (PVL - V_0)/gr
    elseif onset_opt == PVT_model_opt 
        PT = (PVL - VL_LOD_Kissler)/gr
        tp = (1 + onset_frac)*PT
    end

    return PVL, tp, log(10)*gr, log(10)*dr
end