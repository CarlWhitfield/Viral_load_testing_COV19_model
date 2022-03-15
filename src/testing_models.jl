"""
    testing_models.jl

Written by Carl Whitfield, University of Manchester, 2021

Only to be included by definitions.jl
"""

#===================== Testing vs. VL params =====================#

#====> PCR sens (Smith et al. JCMB - fitted) <====#
const PCR_VL0 = 8.522/4.408
const PCR_VLdisp = 4.408
const VL_LOD_PCR = 1.0      #Ct cutoff value for +ve test -- assumed from UKHSA data

#====> LFD params (Porton Down) <====#
const LFD_VLmean = 10.836/2.680
const LFD_VLstd = 1.0/2.680
const LFD_sens_max = 0.75 #assumed

#====> LFD params (Care Home Data) <====#
const VL_lims_2021 = collect(2.0:1.0:7.0)
const LFD_relsens_2021 = [0.004,0.025,0.087,0.162,0.414,0.651,0.8]
const VL_lims_Oct = [4.0,6.0]
const LFD_relsens_Oct = [0.032, 0.255, 0.477]

const PCRTATdist = Gamma(6.75,20/3)

#=================================================================#

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
    if VL > VL_LOD_PCR
        return logistic_function(VL, PCR_VL0, PCR_VLdisp, PCR_sens_max)
    else
        return 0
    end
end

"""
    LFDtest_positive_prob(VL::Float64)

LFD test positive probability given a viral load

## Arguments: 
`VL` = viral load value

## Returns: 
`Float64` = test positive probability
"""
function LFDtest_positive_prob(VL::Float64)
    #test senstivity defined relative to concurrent PCR    
    PCR_prob = PCRtest_positive_prob(VL)
    if LFD_model == SC_data_2021
        VLcat = sum(VL .>= VL_lims_2021) 
        return LFD_relsens_2021[1+VLcat]*PCR_prob
    elseif LFD_model == SC_data_Oct
        VLcat = sum(VL .>= VL_lims_Oct) 
        return LFD_relsens_Oct[1+VLcat]*PCR_prob
    elseif LFD_model == porton_down
        return probit_function(VL, LFD_VLmean, LFD_VLstd, LFD_sens_max)
    end
end

# """ old version
#     LFDtest_positive_prob1(VL::Float64, sens_rel::Float64=1.0)

# LFD test positive probability given a viral load

# ## Arguments: 
# `VL` = viral load value
# `sens_rel` = relative sensitivity scaling factor (default 1)

# ## Returns: 
# `Float64` = test positive probability
# """
# function LFDtest_positive_prob(VL::Float64, sens_rel::Float64=1.0)
#     #sens rel is 1 in Peto et al, could scale this for self-administering effects (sens rel < 1)
#     return sens_rel*probit_function(VL, LFD_VLmean, LFD_VLstd, LFD_sens_max)
# end

function get_pos_profile(sim::Dict, ip::Int, protocol::Int)
    sens_scale = 1.0
    Do_test_prob = 1.0
    if haskey(sim, "P_do_test")
        Do_test_prob = sim["P_do_test"]
    end
    if haskey(sim,"test_sens_scale")
        sens_scale = sim["test_sens_scale"]
    end
    if protocol == PCR_mass_protocol
        sim["test_pos_profiles"][ip] = sens_scale .* Do_test_prob .* PCRtest_positive_prob.(sim["VL_profiles"][ip])
    elseif protocol == LFD_mass_protocol
        sim["test_pos_profiles"][ip] = sens_scale .* Do_test_prob .* LFDtest_positive_prob.(sim["VL_profiles"][ip])
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
function init_testing!(sim::Dict, testing_params::Dict, i_day::Int, Ndays::Int; fill_pos_profiles::Bool=true,
                       different_start::Bool=false)
    #add test positivity profiles to simulation Dict, i_day is start day, Ndays is total length of sim
    #testing params contains "tperiod" (days between test)
    #testing params contains "protocol" (LFD_mass_protocol or PCR_mass_protocol)
    #testing params contains "testing_enforced" (if true all tests are done)
    #test_miss_prob: if "testing_enforced" is false, random probability of missing a test
    #optional: sens_rel: relative scaling of test sensitivity (does not effect spec)
    #optional: policy_adherence: if "testing_enforced" is false, independent probability of adhering to testing
    #(has same effect as sens_rel on test sensitivity, but also affect false positives)
    sim["test_protocol"] = testing_params["protocol"]
    if haskey(testing_params,"sens_rel")
        sim["test_sens_scale"] = testing_params["sens_rel"]
    end
    if (testing_params["testing_enforced"] == false) && haskey(testing_params,"test_miss_prob")
        sim["P_do_test"] = 1.0 - testing_params["test_miss_prob"]
    end
    
    sim["will_isolate_with_test"] = ones(Bool,sim["Ntot"])
    if testing_params["testing_enforced"] == false  #if not enforced, we assume same group who does not isolate with symps does not isolate with tests (simplification)
        if haskey(testing_params,"policy_adherence")
            non_adherers = randsubseq(1:sim["Ntot"], 1 - testing_params["policy_adherence"])
            sim["will_isolate_with_test"][non_adherers] .= false
        else
            sim["will_isolate_with_test"][sim["non_isolators"]] .= false
        end
    end
    sim["testing_paused"] = zeros(Bool,sim["Ntot"])
    sim["resume_testing"] = -ones(Int64,sim["Ntot"])
    
    if different_start
        test_days0 = rand(1:testing_params["tperiod"],sim["Ntot"])
        test_days = Array{Array{Int64,1},1}(undef,sim["Ntot"])
        test_day_counter = ones(Int64,sim["Ntot"])
        for i in 1:sim["Ntot"]
            test_days[i] = collect(test_days0[i]:testing_params["tperiod"]:Int64(ceil(Ndays)))
            test_days[i] = push!(test_days[i], test_days[i][end] + testing_params["tperiod"])
            test_day_counter[i] = 1 + sum(test_days[i] .< i_day)
        end
    else
        test_day0 = rand(1:testing_params["tperiod"])
        test_days = collect(test_day0:testing_params["tperiod"]:Int64(ceil(Ndays)))
        test_days = push!(test_days, test_days[end] + testing_params["tperiod"])
        test_day_counter = 1 + sum(test_days .< i_day)
    end
    sim["test_pos_profiles"] = Array{Array{Float64,1},1}(undef,sim["Ntot"])
    
    if fill_pos_profiles
        sim["test_pos_profiles"] .= get_pos_profile.(Ref(sim), 1:sim["Ntot"], 
                        testing_params["protocol"])
    end
    
    return test_days, test_day_counter
end

function draw_PCR_delays(N::Int)
    #from distribution in social care setting
    delays_hrs = PCR_TaT_scale * rand(PCRTATdist,N)
    delays_days = Int64.(round.(delays_hrs./24))
    #round delay to nearest day
    return delays_days
end

function get_test_probs(VL_profile::Array{Float64,1}, TestDays::Array{Int64,1}, 
                        TestTypes::Array{Int64,1}; Do_PCR_test_prob::Float64 = 1.0,
                        Do_LFD_test_prob::Float64 = 1.0, all_or_nothing::Bool = false)
    
    tp = zeros(length(TestDays))
    TTcondPCR = (TestTypes .== 0) .* (TestDays .<= length(VL_profile)) #PCR is test type 0
    TTcondLFD = (TestTypes .== 1) .* (TestDays .<= length(VL_profile))  #LFD is test type 1
    if all_or_nothing
        Do_all_PCR = (rand() < Do_PCR_test_prob)
        Do_all_LFD = (rand() < Do_LFD_test_prob)
        if Do_all_PCR
            tp[TTcondPCR] = PCRtest_positive_prob.(VL_profile[TestDays[TTcondPCR]])
        else
            tp[TTcondPCR] .= 0
        end
        if Do_all_LFD
            tp[TTcondLFD] = LFDtest_positive_prob.(VL_profile[TestDays[TTcondLFD]])
        else
            tp[TTcondLFD] .= 0 
        end
    else
        tp[TTcondPCR] = Do_PCR_test_prob*PCRtest_positive_prob.(VL_profile[TestDays[TTcondPCR]])
        tp[TTcondLFD] = Do_LFD_test_prob*LFDtest_positive_prob.(VL_profile[TestDays[TTcondLFD]])
    end
    
    return tp
end

function run_testing_scenario!(isol_days::Array{Int64,1}, inf_profile::Array{Float64,1}, 
        test_pos::Array{Float64,1}, test_result_day::Array{Int64,1}, symp_day::Int64, symp_isol::Bool, 
        VL_profile::Array{Float64,1}, Conf_PCR::Array{Bool,1}, 
        Preisolation::Array{Int64,1} = zeros(Int64,0); Day5release::Bool = false, 
        Conf_PCR_fixed::Bool = false, Conf_PCR_sens = PCR_sens_max)
    
    isol_start_day = copy(test_result_day)
    isol_start_prob = copy(test_pos)
    isol_conf_PCR = copy(Conf_PCR)
    if symp_isol && symp_day > 0
        push!(isol_start_day, symp_day)
        push!(isol_start_prob, 1.0)
        push!(isol_conf_PCR, false)   #conf PCR not needed for symptoms
        itd = sortperm(isol_start_day)
        isol_conf_PCR = isol_conf_PCR[itd]
        isol_start_day = isol_start_day[itd]
        isol_start_prob = isol_start_prob[itd]
    end
    go = true
    it = 1
    empty!(isol_days)
    if length(Preisolation) > 0
        push!(isol_days,Preisolation...)
    end
    while go && it <= length(isol_start_day) && isol_start_day[it] <= length(inf_profile)
        if rand() < isol_start_prob[it]  #positive test or symp
            isol_end = isol_start_day[it] + 10
            if Day5release #test daily from day 4, released as soon as 2 consecutive are negative
                DailyLFDResult = zeros(Bool, 6)
                #assume negative once VL profile has ended
                if length(VL_profile) > isol_start_day[it] + 4
                    TestEnd = min(isol_start_day[it] + 9, length(VL_profile))
                    DailyLFDProb = LFDtest_positive_prob.(VL_profile[isol_start_day[it] + 4:TestEnd])
                    DailyLFDResult[1:length(DailyLFDProb)] = rand(length(DailyLFDProb)) .< DailyLFDProb
                end
                EndEarly = false
                k = 5
                while EndEarly == false && k < 10
                    if DailyLFDResult[k - 4] == false && DailyLFDResult[1 + k - 4] == false
                        isol_end = isol_start_day[it] + k
                        EndEarly = true
                    end
                    k = k+1
                end
            end
            if isol_conf_PCR[it]  #assume confirmatory PCR is done immediately
                if isol_start_day[it] <= length(VL_profile)
                    CPCRprob = 0
                    if Conf_PCR_fixed
                        CPCRprob = Conf_PCR_sens
                    else
                        CPCRprob = PCRtest_positive_prob(VL_profile[isol_start_day[it]])  
                    end
                    if rand() > CPCRprob
                        #if conf PCR is negative
                        isol_end = isol_start_day[it] + draw_PCR_delays(1)[1]
                    else
                        go = false   #otherwise isolate in full
                    end
                else
                    go = false   #otherwise isolate in full
                end
            else
                go = false
            end
            if isol_end >= isol_start_day[it]
                push!(isol_days, collect(isol_start_day[it]:1:(isol_end-1))...)
            end
        end
        it = it + 1
    end
    Nisol_days = sum(isol_days)   #record actual number of isolation days
    
    return Nisol_days
end