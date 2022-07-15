include("viral_load_infectivity_testpos.jl")

const scen_names = ["(b) Status Quo","(c1) Fortnightly concurrent PCR","(c2) Fortnightly random PCR", "(d) 3 LFDs per week","(e) 2 LFDs per week","(f) Daily LFDs","(g) Daily LFDs + PCR","(h) 3 LFDs + PCR",
    "(a) No testing"]

#define shift patterns
const shift_pattern = [true, true, true, true, false, false, false, 
                       true, true, true, true, true, false, false]
const twolfd_pattern = [true, false, true, false, false, false, false,
                        true, false, false, true, false, false, false]
const threelfd_pattern = [true, false, true, true, false, false, false,
                          true, false, true, false, true, false, false]

function scenario_1_setup(Ndays::Int, shift_pattern_start::Int)   #2 LFDs + 1 concurrent PCR
    #PCRs happen from first Monday
    #infection day is DayOfWeek[shift_pattern_start]
    #Next Monday is 1 + mod(15 - shift_pattern_start,7)
    #i.e. if shift_pattern_start = 1, infection day is a Monday, so PCR starts on day 1
    #if shift_pattern_start = 2, infection day is a Tuesday, so PCR starts on day 7 etc.
    test_day0 = 1 + mod(15 - shift_pattern_start,7)
    
    Nrepeats = Int64(ceil((Ndays - 13 + shift_pattern_start)/14))
    LFD_test_days_bool = vcat(twolfd_pattern[shift_pattern_start:end], repeat(twolfd_pattern,Nrepeats))
    LFD_test_days_bool = LFD_test_days_bool[1:Ndays]
    
    LFD_test_days = collect(1:length(LFD_test_days_bool))[LFD_test_days_bool]
    PCR_test_days = collect(test_day0:7:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    TestDays = vcat(LFD_test_days, PCR_test_days)
    ResultDays = vcat(LFD_test_days, PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days)), zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end


function scenario_2_setup(Ndays::Int, shift_pattern_start::Int)   #2 LFDs + 1 concurrent fortnightly PCR
    test_day0 = 1 + mod(15 - shift_pattern_start,14)
    
    Nrepeats = Int64(ceil((Ndays - 13 + shift_pattern_start)/14))
    LFD_test_days_bool = vcat(twolfd_pattern[shift_pattern_start:end], repeat(twolfd_pattern,Nrepeats))
    LFD_test_days_bool = LFD_test_days_bool[1:Ndays]
    
    LFD_test_days = collect(1:length(LFD_test_days_bool))[LFD_test_days_bool]
    PCR_test_days = collect(test_day0:14:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    TestDays = vcat(LFD_test_days, PCR_test_days)
    ResultDays = vcat(LFD_test_days, PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days)), zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end

function scenario_3_setup(Ndays::Int, shift_pattern_start::Int)  #2 LFDs + 1 random fortnightly PCR
    Nrepeats = Int64(ceil((Ndays - 13 + shift_pattern_start)/14))
    shifts = vcat(shift_pattern[shift_pattern_start:end], repeat(shift_pattern,Nrepeats))
    shift_days = collect(1:length(shifts))[shifts]
    test_day0 = rand(shift_days[shift_days .< 15]) #random fortnightly PCR on work day
    
    LFD_test_days_bool = vcat(twolfd_pattern[shift_pattern_start:end], repeat(twolfd_pattern,Nrepeats))
    LFD_test_days_bool = LFD_test_days_bool[1:Ndays]
    
    LFD_test_days = collect(1:length(LFD_test_days_bool))[LFD_test_days_bool]
    PCR_test_days = collect(test_day0:14:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    TestDays = vcat(LFD_test_days, PCR_test_days)
    ResultDays = vcat(LFD_test_days, PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days)), zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end

function scenario_4_setup(Ndays::Int, shift_pattern_start::Int)  #3 LFDs
    Nrepeats = Int64(ceil((Ndays - 13 + shift_pattern_start)/14))
    LFD_test_days_bool = vcat(threelfd_pattern[shift_pattern_start:end], repeat(threelfd_pattern,Nrepeats))
    LFD_test_days_bool = LFD_test_days_bool[1:Ndays]
    
    LFD_test_days = collect(1:length(LFD_test_days_bool))[LFD_test_days_bool]
    TestDays = copy(LFD_test_days)
    ResultDays = copy(LFD_test_days)
    TestTypes = ones(length(LFD_test_days))
    itd = sortperm(ResultDays)
    
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end

function scenario_5_setup(Ndays::Int, shift_pattern_start::Int)  #2 LFDs
    Nrepeats = Int64(ceil((Ndays - 13 + shift_pattern_start)/14))
    LFD_test_days_bool = vcat(twolfd_pattern[shift_pattern_start:end], repeat(twolfd_pattern,Nrepeats))
    LFD_test_days_bool = LFD_test_days_bool[1:Ndays]
    
    LFD_test_days = collect(1:length(LFD_test_days_bool))[LFD_test_days_bool]
    TestDays = copy(LFD_test_days)
    TestTypes = ones(length(LFD_test_days))
    
    return TestDays, TestTypes, TestDays
end

function scenario_6_setup(Ndays::Int, shift_pattern_start::Int)  #Daily LFDs
    Nrepeats = Int64(ceil((Ndays - 13 + shift_pattern_start)/14))
    LFD_test_days_bool = vcat(shift_pattern[shift_pattern_start:end], repeat(shift_pattern,Nrepeats))
    LFD_test_days_bool = LFD_test_days_bool[1:Ndays]
    
    TestDays = collect(1:Ndays)[LFD_test_days_bool]
    TestTypes = ones(length(TestDays))
    return TestDays, TestTypes, TestDays #result days same as test days
end

function scenario_7_setup(Ndays::Int, shift_pattern_start::Int)   #Daily LFDs + 1 concurrent PCR
    test_day0 = 1 + mod(15 - shift_pattern_start,7)
    
    Nrepeats = Int64(ceil((Ndays - 13 + shift_pattern_start)/14))
    LFD_test_days_bool = vcat(shift_pattern[shift_pattern_start:end], repeat(shift_pattern,Nrepeats))
    LFD_test_days_bool = LFD_test_days_bool[1:Ndays]
    LFD_test_days = collect(1:Ndays)[LFD_test_days_bool]
    
    PCR_test_days = collect(test_day0:7:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    TestDays = vcat(LFD_test_days, PCR_test_days)
    ResultDays = vcat(LFD_test_days, PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days)), zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end

function scenario_8_setup(Ndays::Int, shift_pattern_start::Int)   #3 LFDs + 1 concurrent PCR
    test_day0 = 1 + mod(15 - shift_pattern_start,7)
    
    Nrepeats = Int64(ceil((Ndays - 13 + shift_pattern_start)/14))
    LFD_test_days_bool = vcat(threelfd_pattern[shift_pattern_start:end], repeat(threelfd_pattern,Nrepeats))
    LFD_test_days_bool = LFD_test_days_bool[1:Ndays]
    
    LFD_test_days = collect(1:length(LFD_test_days_bool))[LFD_test_days_bool]
    PCR_test_days = collect(test_day0:7:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    TestDays = vcat(LFD_test_days, PCR_test_days)
    ResultDays = vcat(LFD_test_days, PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days)), zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end

"""
    init_testing_random!(sim::Dict, testing_params::Dict, i_day::Int, Ndays::Int)

Function to initialise test positive profiles and test isolation probabilities

## Arguments:
`sim` = Dict generated by `init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)`

`testing_params` = Dict containing testing options, must have:

    "scenario" => `String` = Choose from '2_LFD_+_concPCR', '3_LFD', '2_LFD_+_Fort_concPCR'
                                          '2_LFD_+_Fort_randPCR', '2_LFD'

    "comply_prob" => `Float64` = Chance of performing next LFD test (random)

`i_day` = Day to start simulation on (counting from 1)

`Ndays` = Total number of days to simulate (counting from 1)

## Returns:
`Array{Int64,1}` = Days when testing will occur

`Int64` = Next testing day index (from `i_day`)

## See also: 
`init_VL_and_infectiousness(Ntot::Int, Pisol::Float64)`
"""
function init_testing_random!(sim::Dict, testing_params::Dict, Conf_PCR::Bool)
    #add test positivity profiles to simulation Dict, i_day is start day, Ndays is total length of sim
    #testing params contains "tperiod" (days between test)
    #testing params contains "protocol" (LFD_mass_protocol or PCR_mass_protocol)
    #optional: sens_rel: relative sensitivity of LFD as VL -> infinity
    Ndays = length.(sim["VL_profiles"])
    sim["will_isolate_with_test"] = ones(Bool,sim["Ntot"])    
    sim["testing_paused"] = zeros(Bool,sim["Ntot"])
    sim["resume_testing"] = -ones(Int64,sim["Ntot"])
    sim["test_days"] = Array{Array{Int64,1},1}(undef,sim["Ntot"])
    sim["test_types"] = Array{Array{Int64,1},1}(undef,sim["Ntot"])
    sim["test_result_days"] = Array{Array{Int64,1},1}(undef,sim["Ntot"])
    sim["conf_PCR"] = Array{Array{Bool,1},1}(undef,sim["Ntot"])
    TestOutput = []
    if testing_params["scenario"] == scen_names[1]
        TestOutput = scenario_1_setup.(Ndays, sim["shift_pattern_start_day"])
    elseif testing_params["scenario"] == scen_names[2]
        TestOutput = scenario_2_setup.(Ndays, sim["shift_pattern_start_day"])
    elseif testing_params["scenario"] == scen_names[3]
        TestOutput = scenario_3_setup.(Ndays, sim["shift_pattern_start_day"])
    elseif testing_params["scenario"] == scen_names[4]
        TestOutput = scenario_4_setup.(Ndays, sim["shift_pattern_start_day"])
    elseif testing_params["scenario"] == scen_names[5]  
        TestOutput = scenario_5_setup.(Ndays, sim["shift_pattern_start_day"])
    elseif testing_params["scenario"] == scen_names[6]  
        TestOutput = scenario_6_setup.(Ndays, sim["shift_pattern_start_day"])   
    elseif testing_params["scenario"] == scen_names[7]  
        TestOutput = scenario_7_setup.(Ndays, sim["shift_pattern_start_day"])    
    elseif testing_params["scenario"] == scen_names[8]  
        TestOutput = scenario_8_setup.(Ndays, sim["shift_pattern_start_day"])       
    end
    for i in 1:sim["Ntot"]
        sim["test_days"][i] = TestOutput[i][1]
        sim["test_types"][i] = TestOutput[i][2]
        sim["test_result_days"][i] = TestOutput[i][3]
        sim["conf_PCR"][i] = zeros(Bool,length(sim["test_days"][i]))
        if Conf_PCR
            sim["conf_PCR"][i][sim["test_types"][i] .== 1] .= true
        end
    end
    
    sim["test_pos_prob"] = get_test_probs.(sim["VL_profiles"], sim["test_days"], 
        sim["test_types"]; Do_LFD_test_prob = testing_params["comply_prob"],
        all_or_nothing = testing_params["aon_compliance"]) 
    #assuming random compliance
end

function get_no_testing_scenario(Ntot::Int, Pisol::Float64; 
        Day5release_bool::Bool = false)
    sim_baseline = init_VL_and_infectiousness(Ntot, Pisol)
    sim_baseline["shift_pattern_start_day"] = rand(1:length(shift_pattern),Ntot)
    sim_baseline["test_pos_prob"] = Array{Array{Float64,1},1}(undef,Ntot)
    sim_baseline["test_result_days"]  = Array{Array{Int64,1},1}(undef,Ntot)
    Conf_PCR_h = fill(zeros(Bool,0),sim_baseline["Ntot"])
    sim_baseline["isol_days"] = Array{Array{Int64,1},1}(undef,Ntot)
    sim_baseline["tot_isol_days"] = zeros(Int64,Ntot)
    for j in 1:Ntot
        sim_baseline["test_pos_prob"][j] = zeros(0)
        sim_baseline["test_result_days"][j] = zeros(Int64,0)
        sim_baseline["isol_days"][j] = zeros(Int64,0)
        sim_baseline["tot_isol_days"][j] = run_testing_scenario!(sim_baseline["isol_days"][j],
            sim_baseline["infection_profiles"][j], sim_baseline["test_pos_prob"][j],
            sim_baseline["test_result_days"][j], sim_baseline["symp_day"][j] .+ 1, 
            sim_baseline["will_isolate"][j], sim_baseline["VL_profiles"][j], 
            Conf_PCR_h[j]; Day5release=Day5release_bool)
    end
    
    return sim_baseline
end

function run_testing_scenarios_vs_baseline(sim_baseline::Dict, LFD_comply::Float64, 
                                      Conf_PCR::Bool; Day5release_bool::Bool = false,
                                      LFD_AllorNone::Bool = false)
    Ntot = sim_baseline["Ntot"]
    Nscens = length(scen_names)
    sim_scens = Array{Dict,1}(undef,Nscens)
    for i in 1:(Nscens-1)
        sim_scens[i] = copy(sim_baseline)
        init_testing_random!(sim_scens[i], Dict("scenario"=>scen_names[i],
                 "comply_prob"=>LFD_comply, "aon_compliance"=>LFD_AllorNone), Conf_PCR)
        
        sim_scens[i]["isol_days"] = copy.(sim_baseline["isol_days"])
        sim_scens[i]["tot_isol_days"] = run_testing_scenario!.(sim_scens[i]["isol_days"],
            sim_scens[i]["infection_profiles"],  sim_scens[i]["test_pos_prob"], sim_scens[i]["test_result_days"], 
            sim_scens[i]["symp_day"] .+ 1,  sim_scens[i]["will_isolate"], sim_scens[i]["VL_profiles"], 
            sim_scens[i]["conf_PCR"]; Day5release=Day5release_bool)
    end
    sim_scens[Nscens] = copy(sim_baseline)
    for i in 1:Nscens
        sim_scens[i]["inf_days"] = zeros(Int64, Ntot)
        sim_scens[i]["tot_work_isol_days"] = zeros(Int64, Ntot)
        sim_scens[i]["inf_profile_isolation"] = copy(sim_scens[i]["infection_profiles"])
        Ndays = length.(sim_scens[i]["VL_profiles"])
        Nrepeats = Int64.(ceil.((Ndays .- 13 .+ sim_scens[i]["shift_pattern_start_day"]) ./ 14))
        for j in 1:Ntot
            if length(sim_scens[i]["isol_days"][j]) > 0
                maxisolday = max(sim_scens[i]["isol_days"][j]...)
                if maxisolday > Ndays[j]
                    Nrepeats[j] = Int64(ceil((maxisolday - 13 + sim_scens[i]["shift_pattern_start_day"][j]) / 14))
                end
            end
            shifts = vcat(shift_pattern[sim_scens[i]["shift_pattern_start_day"][j]:end], 
                          repeat(shift_pattern,Nrepeats[j]))
            inf_shifts = shifts[1:Ndays[j]]
            #count people as not infectious if not at work
            sim_scens[i]["infection_profiles"][j][inf_shifts .== false] .= 0.0
            #extract isolation days that occur in infectious window
            inf_isol_days = sim_scens[i]["isol_days"][j][sim_scens[i]["isol_days"][j] .<= 
                                            length(sim_scens[i]["infection_profiles"][j])]
            #copy over infectious profile with shifts removed
            sim_scens[i]["inf_profile_isolation"][j] = copy(sim_scens[i]["infection_profiles"][j])
            #set isolation days to zero in modified inf profile
            sim_scens[i]["inf_profile_isolation"][j][inf_isol_days] .= 0.0
            sim_scens[i]["inf_days"][j] = sum((sim_scens[i]["inf_profile_isolation"][j] .> 0))
            sim_scens[i]["tot_work_isol_days"][j] = sum(shifts[sim_scens[i]["isol_days"][j]])
        end
    end
    
    return sim_scens, scen_names
end








