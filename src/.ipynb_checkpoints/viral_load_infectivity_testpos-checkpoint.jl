"""
    viral_load_infectivity_testpos.jl

Written by Carl Whitfield, University of Manchester, 2021
"""
# import Pkg
# Pkg.add("DataFrames")
# Pkg.add("CSV")
# Pkg.add("Distributions")
# Pkg.add("Random")
# Pkg.add("StatsBase")
# Pkg.add("SpecialFunctions")
# Pkg.add("Plots")

using DataFrames
using CSV
using Distributions
using Random
using StatsBase
using SpecialFunctions
using Plots

#testing protocol definitions
const PCR_mass_protocol = "PCR_mass_testing"
const LFD_mass_protocol = "LFD_mass_testing"

#VL trajectory params (Kissler)
const beta_onset = 0.7
const alpha_onset = 1.5*beta_onset
const beta_peak = 0.7
const alpha_peak = 3.5*beta_peak
const beta_decay = 0.3
const alpha_decay_symp = 10.5*beta_decay
const alpha_decay_asymp = 6.7*beta_decay
const peakVL_mean = 7.533
const peakVL_sd = 1.164
const VL_LOD_Kissler = 2.66  #starting VL in log10 copies/ml
#PCR swab params (Smith et al. JCMB)
const PCR_VL0 = 8.522/4.408
const PCR_VLdisp = 4.408
const PCR_sens_max = 0.85 #assumed
#LFD params (Porton Down)
const LFD_VLmean = 10.836/2.680
const LFD_VLstd = 1.0/2.680
const LFD_sens_max = 0.75 #assumed
#VL ~ infectivity relation (Marks)
const inf_dep = 1.3
const VL_ref = peakVL_mean
const PIsigma = 0.5
const PImu = -0.5*PIsigma^2
#Viral load where people stop being infectious
const inf_VL_cutoff = 6.0 #3.0 -- should only have minor effect
#this scaling means that, on average, p_inf is approx 0.06 for hour long face-to-face contacts within 14 days of infection.
const j0scale = 4.1 * 1.1 *(4.0/(exp(log(inf_dep)^2*peakVL_sd^2/2)*(1 +
                erf((peakVL_mean - inf_VL_cutoff + log(inf_dep)*peakVL_sd^2)/(sqrt(2)*peakVL_sd)))))
#correlation gradient between peak_inf and pasymp
# const pasymp_vl7 = 0.62
# const pasymp_vl89 = 0.5
# const pasymp_vl10 = 0.34
const p_asymp = 0.5
#const asymp_frac = cdf(Normal(peakVL_mean,peakVL_sd),7.0)*pasymp_vl7 +
#    (cdf(Normal(peakVL_mean,peakVL_sd),9.0) - cdf(Normal(peakVL_mean,peakVL_sd),7.0))*pasymp_vl89 +
#    (1 - cdf(Normal(peakVL_mean,peakVL_sd),9.0))*pasymp_vl10
#relative infectivity of asymptomatics
#symp time distribution
const symp_beta = 4.84 / 2.6^2
const symp_alpha = 4.84 * symp_beta

const PVLdist = Normal(peakVL_mean,peakVL_sd)
const OTdist = Gamma(alpha_onset,1/beta_onset)
const PTdist = Gamma(alpha_peak,1/beta_peak)
const DTsympdist = Gamma(alpha_decay_symp,1/beta_decay)
const DTasympdist = Gamma(alpha_decay_asymp,1.0/beta_decay)

""" 
    generate_peak_viral_load()

Generate and return a single peak viral load value. 

## Returns:
`Float` = peak viral load (log10 copies/ml)
"""
function generate_peak_viral_load()
    return rand(PVLdist)
end

""" 
    generate_onset_time()

Generate and return a single onset time (time from infection
when viral load is detectable).

## Returns: 
`Float` = onset time (days)
"""
function generate_onset_time()
    return rand(OTdist)
end

""" 
    generate_peak_time()

Generate and return a single peak time (time from onset 
to peak viral load).

## Returns: 
`Float` = peak time (days)
"""
function generate_peak_time()
    return rand(PTdist)
end

"""
    generate_asymptomatic()

Generate and return a boolean value determining whether
somebody will develop symptoms.

## Returns: 
`Bool` = whether or not the person in asymptomatic (true = 
          asymptomatic, false = symptoms)
"""
function generate_asymptomatic()
    #randomly generate if people are asymptomatic
    return (rand() < p_asymp)
end

"""
    generate_decay_time(Asymp::Bool)

Generate and return a single decay time (time from peak viral load 
to undetectable VL).

## Arguments: 
`Asymp` = whether or not the person in asymptomatic

## Returns: 
`Float` = Decay time (days)

## See also: 
`generate_asymptomatic()`
"""
function generate_decay_time(Asymp::Bool)
    DT = 
    if Asymp
        DT = rand(DTasympdist)
    else
        DT = rand(DTsympdist)
    end
    return DT
end

"""
    generate_peak_viral_loads(Ntot::Int)

Generate and return multiple peak viral load values.

## Arguments:
`Ntot` = Number of values to generate

## Returns: 
`Array{Float64,1}` = Peak viral loads (log10 copies/ml)

## See also:
`generate_peak_viral_load()`
"""
function generate_peak_viral_loads(Ntot::Int)
    return rand(PVLdist,Ntot)
end

"""
    generate_onset_times(Ntot::Int)

Generate and return multiple onset times (time from infection 
when viral load is detectable).

## Arguments: 
`Ntot` = Number of values to generate

## Returns: 
`Array{Float64,1}` = onset times (days)

## See also:
`generate_onset_time()`
"""
function generate_onset_times(Ntot::Int)
    return rand(OTdist,Ntot)
end

"""
    generate_peak_times(Ntot::Int)

Generate and return multiple onset times (time from infection 
when viral load is detectable).

## Arguments: 
`Ntot` = Number of values to generate

## Returns: 
`Array{Float64,1}` = peak times (days)

## See also:
`generate_peak_time()`
"""
function generate_peak_times(Ntot::Int)
    return rand(PTdist,Ntot)
end

"""
    generate_asymptomatics(Ntot::Int)

Generate and return multiple asymptomatic indicators

## Arguments: 
`Ntot` = Number of values to generate

## Returns: 
`Array{Bool,1}` = Asymptomatic statuses

## See also:
`generate_asymptomatic()`
"""
function generate_asymptomatics(Ntot::Int)
    #randomly generate if people are asymptomatic
    return (rand(Ntot) .< p_asymp)
end

"""
    generate_decay_times(Ntot::Int, Asymp::Array{Bool,1})

Generate and return multiple decay times (time from peak viral 
load to undetectable VL).

## Arguments: 
`Ntot` = Number of values to generate

## Returns: 
`Array{Float64,1}` = Array of decay times

## See also:
`generate_decay_time()`
"""
function generate_decay_times(Ntot::Int, Asymp::Array{Bool,1})
    DTs = generate_decay_time.(Asymp)
    return DTs
end

"""
    generate_peak_inf(peak_VL::Float64)

Generate and return peak infectivity

## Arguments: 
`peak_VL` = Peak viral load (log10 copies/ml)

## Returns: 
`Float64` = Peak infectivity (normalised units)

## See also:
`generate_peak_viral_load()`
"""
function generate_peak_inf(peak_VL::Float64)
    #assume reference person has average infectivity of 1 (so peak of 2)
    value = j0scale * inf_dep^(peak_VL - VL_ref)
    r = rand(LogNormal(PImu,PIsigma))
    return r*value
end

"""
    generate_symp_time(OT::Float64, PT::Float64, DT::Float64)

Generate and return single time of symptom onset, 
bounded by the times where viral load is at detectable levels

## Arguments: 
`OT` = Onset time (days)

`PT` = Time from onset to peak VL (days)

`DT` = Time from peak to undetectable VL (days)

## Returns: 
`Float64` = Symptom onset time (days)

## See also:
`generate_onset_time()`, `generate_peak_time()`, `generate_decay_time()`
"""
function generate_symp_time(OT::Float64,PT::Float64,DT::Float64)
    Strunc = truncated(Gamma(symp_alpha,1.0/symp_beta), OT, OT + PT + DT)
    return rand(Strunc)
end

"""
    build_viral_load_distribution!(v::Array{Float64,1}, inf::Array{Float64,1})

Generate a viral load and infectivity trajectory

## Arguments: 
`sim` = Dict container with simulation info

` i` = Index of individual to build

## See also:
`build_viral_load_distributions!(sim::Dict)`
"""
function build_viral_load_distribution!(sim::Dict, index::Int64)
    PVL = generate_peak_viral_load()
    OT = generate_onset_time()
    PT = generate_peak_time()
    Asymp = generate_asymptomatic()
    DT = generate_decay_time(Asymp)
    PInf = generate_peak_inf(PVL)
    #generate blank containers for viral load and infectivity
    
    v = zeros(Int64(ceil(OT + PT + DT + 10)))
    inf = zeros(Int64(ceil(OT + PT + DT + 10)))
    
    m_up = (PVL - VL_LOD_Kissler)/PT     #slope up to peak VL
    m_down = (PVL - VL_LOD_Kissler)/DT   #slope down from peak VL
    T = length(v)
    i = 1:T                              #day index
    t = i .- 1                           #days since infection
    i1 = Int64(ceil(OT))                 #onset day index
    i2 = Int64(ceil(OT + PT))            #peak VL day index
    i3 = Int64(ceil(OT + PT + DT))       #end VL day index
    cond1 = (i .<= i2)
    #set VL before peak
    v[cond1] = VL_LOD_Kissler .+ m_up .* (t[cond1] .- OT)
    cond2 = (i .> i2)
    #set VL after peak
    v[cond2] = PVL .- m_down .* (t[cond2] .- OT .- PT)
    #generate symptom onset time
    ST = generate_symp_time(OT,PT,DT)
    if PVL > inf_VL_cutoff    #if peak viral load high enough to be infectious
        #infectiousness goes linearly with VL above cutoff
        tinf_start = OT + (inf_VL_cutoff - VL_LOD_Kissler)/m_up
        tinf_end = OT + PT + (PVL - inf_VL_cutoff)/m_down
        cum_inf = zeros(T+1)
        t_inf = collect(-1:(T-1)) .+ 0.5
           
        #ensure cumulative infectiousness is preserved in interpolation
        cond1 = (t_inf .>= tinf_start) .* (t_inf .<  OT + PT)
        cum_inf[cond1] = 0.5 * PInf .* (t_inf[cond1] .- tinf_start).^2 / 
                               (OT + PT - tinf_start)
        
        cond2 = (t_inf .>= OT + PT) .* (t_inf .<=  tinf_end)
        c_inf_peak = 0.5 * PInf * (OT + PT - tinf_start)
        cum_inf[cond2] = c_inf_peak .+ PInf .* (t_inf[cond2] .- OT .- PT) .*
                    (1.0 .- 0.5 .* (t_inf[cond2] .- OT .- PT) ./ (tinf_end  - OT - PT))
        cum_inf[t_inf .>  tinf_end] .= 0.5*PInf*(tinf_end - tinf_start)
        
        inf .= cum_inf[2:(T+1)] .- cum_inf[1:T]
        nel = 1:length(t)
        deleteat!(inf, nel[t .> Int64(ceil(tinf_end))])
    end
    SD = Int64(round(rand() + ST))

    sim["symp_day"][index] = SD
    sim["VL_mag"][index] = PVL
    sim["asymptomatic"][index] = Asymp
    sim["VL_profiles"][index] = v
    sim["infection_profiles"][index] = inf
end

"""
    build_viral_load_distributions!(sim::Dict)

Generate multiple viral load and infectivity trajectories.

## Arguments:
`sim` = Dict container with the following keys:

    "VL_profiles"=>Array{Array{Float64,1},1} = Each entry i is
        a container for a generated viral load trajectory for
        individual i, interpolated to one value per day.

    "infection_profiles"=>Array{Array{Float64,1},1} = Each 
        entry i is a container for a generated infectivity profile
        of individual i, interpolated to one value per day.

    "symp_day"=>Array{Int64,1} = Each entry i is the symptom
        onset time of individual i.

    "VL_mag"=>Array{Float64,1} = Each entry i is the peak viral
        load of individual i.

    "asymptomatic"=>Array{Float64,1}  = Each entry i is the 
        asymptomatic status of individual i.

    All arrays must have the same size (number of trajectories
    to be generated), and the entries will be overwritten with
    the generated values.

## Returns: 
none

## See also
`init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)`
"""
function build_viral_load_distributions!(sim::Dict)
    build_viral_load_distribution!.(Ref(sim), collect(1:sim["Ntot"]))
end

# function infectivity(VL::Array{Float64,1}, peak_inf::Float64, peak_VL::Float64)
#     j = zeros(length(VL))
#     t = 0:(length(VL)-1)
#     cond1 = (VL .> inf_VL_cutoff)
#     j[cond1] = peak_inf .* (VL[cond1] .-  inf_VL_cutoff) ./ (peak_VL - inf_VL_cutoff)
#     return j
# end

# # function infectivity_alt(VL::Array{Float64,1}, peak_inf::Float64, peak_VL::Float64, peak_VL_time::Float64, decay_time::Float64)
# #     #assume infectivity scales linearly over time with log10 copies/ml above cutoff
# #     j = zeros(length(VL))
# #     t = 0:(length(VL)-1)
# #     cond1 = (VL .> inf_VL_cutoff)
# #     j[cond1] = peak_inf .* (VL[cond1] .- inf_VL_cutoff) ./ (peak_VL - inf_VL_cutoff)
# #     cond2 = (VL .> inf_VL_cutoff) .* (t .>= peak_VL_time)
# #     j[cond2] = peak_inf .* (1 .+ (peak_VL_time .- t[cond2])) ./ (0.38 * decay_time)

# #     return j
# # end

"""
    logistic_function(x::Float64, x0::Float64, k::Float64, ymax::Float64)

Scaled logistic function 

## Arguments: 
`x`

`x0` = Point where logistic function is at 50% of max value

`k` = Slope of logistic function

`ymax` = Maximum value of logistic function

## Returns: 
`Float64` = ymax / (1  + exp(-k*(x-x0)))
"""
function logistic_function(x::Float64, x0::Float64, k::Float64, ymax::Float64)
    return (ymax / (1  + exp(-k*(x-x0))))
end

"""
    probit_function(x::Float64, x0::Float64, xstd::Float64, ymax::Float64)

Scaled probit function 

## Arguments: 
`x`

`x0` = Mean of normal distribution kernel

`xstd` = Standar deviation of normal distribution kernel

`ymax` = Maximum value of probit function

## Returns: 
`Float64` = ymax * cdf(Normal(x0,xstd),x)
"""
function probit_function(x::Float64, x0::Float64, xstd::Float64, ymax::Float64)
    return ymax * cdf(Normal(x0,xstd),x)
end

"""
    PCRtest_positive_prob(VL::Float64)

PCR test positive probability given a viral load

## Arguments: 
`VL` = viral load value (log10 copies/ml)

## Returns: 
`Float64` = test positive probability
"""
function PCRtest_positive_prob(VL::Float64)
    #assume PCR is 100% accurate at detecting viral material on swab below Ct = 40
    if VL > VL_LOD_Kissler
        return logistic_function(VL, PCR_VL0, PCR_VLdisp, PCR_sens_max)
    else
        return 0
    end
end

"""
    LFDtest_positive_prob(VL::Float64, sens_rel::Float64=1.0)

LFD test positive probability given a viral load

## Arguments: 
`VL` = viral load value
`sens_rel` = relative sensitivity scaling factor (default 1)

## Returns: 
`Float64` = test positive probability
"""
function LFDtest_positive_prob(VL::Float64, sens_rel::Float64=1.0)
    #sens rel is 1 in Peto et al, could scale this for self-administering effects (sens rel < 1)
    return sens_rel*probit_function(VL, LFD_VLmean, LFD_VLstd, LFD_sens_max)
end

"""
    generate_isolations(Ntot::Int, Pisol::Float64)

Randomly generate whether people obey symptomatic isolation or not

## Arguments: 
`Ntot` = Number to generate (indices are 1:Ntot)
`Pisol` = Probability of isolation

## Returns: 
`Array{Int,1}` = Array of indices of those who will isolate at symptom
                 onset (if symptomatic)
"""
function generate_isolations(Ntot::Int, Pisol::Float64)
    return randsubseq(1:Ntot, Pisol)
end

"""
    init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)

Randomly generate all viral load and infectivity data for a simulation

## Arguments:
`Ntot` = Number of individuals to simulate
`Pisol` = Probability of self-isolation at symptom onset

## Returns: 
`Dict` = Container with generated parameters including:

        "Ntot" => `Int` = Number of individuals in simulation (indexed 1:Ntot)

        "asymptomatic" => `Array{Bool,1}` = Asymptomatic status for each individual

        "isolation_time" => `Array{Int64,1}` = Day of isolation for each individual 
                                           (to be filled later)

         "will_isolate" => `Array{Bool,1}` = Indicator of whether an individual will
                                         isolate at symptom onset time
        
         "VL_mag" => `Array{Float64,1}` = Peak viral load for each individual

         "inf_mag" => `Array{Float64,1}` = Peak infectivity for each individual

         "infection_profiles" => `Array{Array{Float64,1},1}` = infectivity for each day
                            since infection, for each individual

         "VL_profiles" => `Array{Array{Float64,1},1}` = viral load for each day
                            since infection, for each individual

         "days_infectious" => `Array{Int64,1}` = Number of days spent with infectivity
                                               > 0 for each individual

          "symp_day" => `Array{Int64,1}` = Day at which individuals would isolate due to
                                         symptom onset

          "non_isolators" => `Array{Int64,1}` = Array of indices ofindividuals who would 
                                          refuse to isolate if they developed symptoms
"""
function init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)
    #simulate Ntot people with baseline isolation probability Pisol
    sim = Dict("Ntot"=>Ntot, "asymptomatic"=>zeros(Bool, Ntot),
               "VL_mag"=>zeros(Float64,Ntot),
               "isolation_time"=>zeros(Int64, Ntot),
               "will_isolate"=>zeros(Bool, Ntot),
               "inf_mag"=>zeros(Float64, Ntot),
               "infection_profiles"=>Array{Array{Float64,1},1}(undef,0),
               "VL_profiles"=>Array{Array{Float64,1},1}(undef,0),
               "days_infectious" => zeros(Int64,Ntot))
    for i in 1:Ntot
        push!(sim["infection_profiles"],zeros(0))
        push!(sim["VL_profiles"],zeros(0))
    end
    build_viral_load_distributions!(sim)
    sim["will_isolate"][generate_isolations(Ntot, Pisol)] .= true
    nr = 1:sim["Ntot"]
    sim["non_isolators"] = nr[(sim["will_isolate"] .== false)]
    sim["will_isolate"][sim["asymptomatic"]] .= false #asymptomatics don't self isolate, but are not "non isolators"

    return sim
end

function get_pos_profile(sim::Dict, ip::Int, protocol::Int; LFD_sens_rel::Float64 = 1.0)
    if protocol == PCR_mass_protocol
        sim["test_pos_profiles"][ip] = PCRtest_positive_prob.(sim["VL_profiles"][ip])
    elseif protocol == LFD_mass_protocol
        sim["test_pos_profiles"][ip] = LFDtest_positive_prob.(sim["VL_profiles"][ip], Ref(sens_rel))
    end
end
    
    
"""
    init_testing!(sim::Dict, testing_params::Dict, i_day::Int, Ndays::Int)

Function to initialise test positive profiles and test isolation probabilities

## Arguments:
`sim` = Dict generated by `init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)`

`testing_params` = Dict containing testing options, must have:

    "new_comply_prob" => `Float` = Probability that people who would not comply
                       with symptom onset isolation, will comply with testing isolation

    "tperiod" => `Int` = Days between repeat tests

    "protocol" => `String` = Current available options are given by the parameters
                            `PCR_mass_protocol` and `LFD_mass_protocol`

`i_day` = Day to start simulation on (counting from 1)

`Ndays` = Total number of days to simulate (counting from 1)

## Returns:
`Array{Int64,1}` = Days when testing will occur

`Int64` = Next testing day index (from `i_day`)

## See also: 
`init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)`
"""
function init_testing!(sim::Dict, testing_params::Dict, i_day::Int, Ndays::Int; fill_pos_profiles::Bool=true)
    #add test positivity profiles to simulation Dict, i_day is start day, Ndays is total length of sim
    #testing params contains "tperiod" (days between test)
    #testing params contains "protocol" (LFD_mass_protocol or PCR_mass_protocol)
    #optional: sens_rel: relative sensitivity of LFD as VL -> infinity
    sim["test_protocol"] = testing_params["protocol"]
    sim["will_isolate_with_test"] = ones(Bool,sim["Ntot"])
    sim["will_isolate_with_test"][sim["non_isolators"]] .= false
    new_compliers = randsubseq(sim["non_isolators"], testing_params["new_comply_prob"])
    if length(new_compliers) > 0
        sim["will_isolate_with_test"][new_compliers] .= true
    end
    test_day0 = rand(1:testing_params["tperiod"])
    test_days = collect(test_day0:testing_params["tperiod"]:Int64(ceil(Ndays)))
    test_days = push!(test_days, test_days[end] + testing_params["tperiod"])
    test_day_counter = 1 + sum(test_days .< i_day)
    sim["test_pos_profiles"] = Array{Array{Float64,1},1}(undef,sim["Ntot"])
    if fill_pos_profiles
        if haskey(testing_params,"sens_rel")
            sim["test_pos_profiles"] .= get_pos_profile.(Ref(sim), 1:sim["Ntot"], 
                        testing_params["protocol"], testing_params["sens_rel"])
        else
            sim["test_pos_profiles"] .= get_pos_profile.(Ref(sim), 1:sim["Ntot"], 
                        testing_params["protocol"])
        end
    end

    return test_days, test_day_counter
end

"""
    plot_infection_profiles(sim::Dict)

Create a plot of all infectiousness profiles

## Arguments:
`sim` = Dict created by `init_VL_and_infectiousness`

## Returns:

## See also: 
`init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)`
"""
function plot_infection_profiles(sim::Dict)
    if sim["Ntot"] > 0
        pp = Plots.plot(sim["infection_profiles"][1], xlabel = "Days", ylabel = "Infectiousness")
        for i in 2:sim["Ntot"]
            pp = Plots.plot!(sim["infection_profiles"][i])
        end
    end

    pout = Plots.display(pp)
    
    return pout
end