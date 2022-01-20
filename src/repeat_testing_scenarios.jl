include("viral_load_infectivity_testpos_v2.jl")

const scen_names = ["(b) Status Quo","(c1) Fortnightly concurrent PCR","(c2) Fortnightly random PCR", "(d) 3 LFDs per week","(e) 2 LFDs per week","(f) Daily LFDs","(g) Daily LFDs + PCR","(h) 3 LFDs + PCR",
    "(a) No testing"]

use_ke_model = false

function scenario_1_setup(Ndays::Int)   #2 LFDs + 1 concurrent PCR
    #assume day0 is random for each person
    test_day0 = rand(1:7)
    LFD_test_days1 = collect((test_day0+3):7:Ndays)
    LFD_test_days2 = collect(test_day0:7:Ndays)    
    PCR_test_days = collect(test_day0:7:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    if test_day0 > 4
        push!(LFD_test_days1,test_day0-4)
    end
    TestDays = vcat(LFD_test_days1,LFD_test_days2,PCR_test_days)
    ResultDays = vcat(LFD_test_days1,LFD_test_days2,PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days1) + length(LFD_test_days2)),
                             zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end


function scenario_2_setup(Ndays::Int)   #2 LFDs + 1 concurrent fortnightly PCR
    #assume day0 is random for each person
    test_day0 = rand(1:14)
    LFD_test_days1 = collect((test_day0+3):7:Ndays)
    LFD_test_days2 = collect(test_day0:7:Ndays)     
    PCR_test_days = collect(test_day0:14:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    if test_day0 > 4
        push!(LFD_test_days1,test_day0-4)
        if test_day0 > 7
            push!(LFD_test_days1,test_day0-7)
            if test_day0 > 11
                push!(LFD_test_days1,test_day0-11)
            end
        end
    end
    TestDays = vcat(LFD_test_days1,LFD_test_days2,PCR_test_days)
    ResultDays = vcat(LFD_test_days1,LFD_test_days2,PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days1) + length(LFD_test_days2)),
                             zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end

function scenario_3_setup(Ndays::Int)  #2 LFDs + 1 random fortnightly PCR
    LFDtest_day0 = rand(1:7)
    PCRtest_day0 = rand(1:14)
    LFD_test_days1 = collect((LFDtest_day0+3):7:Ndays)
    LFD_test_days2 = collect(LFDtest_day0:7:Ndays)        
    PCR_test_days = collect(PCRtest_day0:14:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    if LFDtest_day0 > 4
        push!(LFD_test_days1,LFDtest_day0-4)
    end 
    TestDays = vcat(LFD_test_days1,LFD_test_days2,PCR_test_days)
    ResultDays = vcat(LFD_test_days1,LFD_test_days2,PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days1) + length(LFD_test_days2)),
                             zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end

function scenario_4_setup(Ndays::Int)  #3 LFDs
    test_day0 = rand(1:7)
    LFD_test_days1 = collect((test_day0+4):7:Ndays)
    LFD_test_days2 = collect((test_day0+2):7:Ndays) 
    LFD_test_days3 = collect(test_day0:7:Ndays)
    if test_day0 > 3
        push!(LFD_test_days2,test_day0-3)
        if test_day0 > 5
            push!(LFD_test_days1,test_day0-5)
        end
    end 
    TestDays = vcat(LFD_test_days1,LFD_test_days2,LFD_test_days3)
    TestTypes = ones(length(TestDays))
    itd = sortperm(TestDays)
    return TestDays[itd], TestTypes[itd], TestDays[itd] #result days same as test days
end

function scenario_5_setup(Ndays::Int)  #2 LFDs
    test_day0 = rand(1:7)
    LFD_test_days1 = collect((test_day0+3):7:Ndays)
    LFD_test_days2 = collect(test_day0:7:Ndays)
    if test_day0 > 4
        push!(LFD_test_days1,test_day0-4)
    end 
    TestDays = vcat(LFD_test_days1,LFD_test_days2)
    TestTypes = ones(length(TestDays))
    itd = sortperm(TestDays)
    return TestDays[itd], TestTypes[itd], TestDays[itd] #result days same as test days
end

function scenario_6_setup(Ndays::Int)  #7 LFDs
    TestDays = collect(1:Ndays)
    TestTypes = ones(length(TestDays))
    return TestDays, TestTypes, TestDays #result days same as test days
end

function scenario_7_setup(Ndays::Int)   #Daily LFDs + 1 concurrent PCR
    #assume day0 is random for each person
    test_day0 = rand(1:7)
    LFD_test_days = collect(1:Ndays)  
    PCR_test_days = collect(test_day0:7:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    TestDays = vcat(LFD_test_days,PCR_test_days)
    ResultDays = vcat(LFD_test_days,PCR_result_days)
    TestTypes = vcat(ones(length(LFD_test_days)), zeros(length(PCR_test_days)))
    itd = sortperm(ResultDays)
    return TestDays[itd], TestTypes[itd], ResultDays[itd]
end

function scenario_8_setup(Ndays::Int)   #3 LFDs + 1 concurrent PCR
    #assume day0 is random for each person
    test_day0 = rand(1:7)
    
    LFD_test_days1 = collect((test_day0+4):7:Ndays)
    LFD_test_days2 = collect((test_day0+2):7:Ndays) 
    LFD_test_days3 = collect(test_day0:7:Ndays)
    if test_day0 > 3
        push!(LFD_test_days2,test_day0-3)
        if test_day0 > 5
            push!(LFD_test_days1,test_day0-5)
        end
    end 
    LFD_test_days = vcat(LFD_test_days1,LFD_test_days2,LFD_test_days3)
    
    PCR_test_days = collect(test_day0:7:Ndays)
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    
    TestDays = vcat(LFD_test_days,PCR_test_days)
    ResultDays = vcat(LFD_test_days,PCR_result_days)
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
        TestOutput = scenario_1_setup.(Ndays)
    elseif testing_params["scenario"] == scen_names[2]
        TestOutput = scenario_2_setup.(Ndays)
    elseif testing_params["scenario"] == scen_names[3]
        TestOutput = scenario_3_setup.(Ndays)
    elseif testing_params["scenario"] == scen_names[4]
        TestOutput = scenario_4_setup.(Ndays)
    elseif testing_params["scenario"] == scen_names[5]  
        TestOutput = scenario_5_setup.(Ndays)
    elseif testing_params["scenario"] == scen_names[6]  
        TestOutput = scenario_6_setup.(Ndays)   
    elseif testing_params["scenario"] == scen_names[7]  
        TestOutput = scenario_7_setup.(Ndays)    
    elseif testing_params["scenario"] == scen_names[8]  
        TestOutput = scenario_8_setup.(Ndays)       
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

"""
    run_testing_scenarios_impact(Ntot::Int, Pisol::Float64, LFD_comply::Float64, 
                                      Conf_PCR::Bool; )

## Arguments: 
`Ntot` = Number of profiles to generate
`LFD_comply` = probability of doing each LFD test
`Conf_PCR` = Do confirmatory PCR after positive LFD (if negative stop isolating early)

## Returns: 
`sim_scens` = Array of Dicts, one for each scenario simulated
`sim_names` = Names of scenarios simulated
 """
function run_testing_scenarios_impact(Ntot::Int, Pisol::Float64, LFD_comply::Float64, 
                            Conf_PCR::Bool; Day7release_bool::Bool = false,
                            Day67tests_bool::Bool = true, LFD_AllorNone::Bool = false)
    sim_baseline = init_VL_and_infectiousness(Ntot, Pisol)
    Nscens = length(scen_names)
    sim_scens = Array{Dict,1}(undef,Nscens)
    for i in 1:(Nscens-1)
        sim_scens[i] = copy(sim_baseline)
        init_testing_random!(sim_scens[i], Dict("scenario"=>scen_names[i],
                    "comply_prob"=>LFD_comply, "aon_compliance"=>LFD_AllorNone), Conf_PCR)
        
        sim_scens[i]["inf_profile_isolation"] = copy.(sim_scens[i]["infection_profiles"])
        sim_scens[i]["isol_days"] = run_testing_scenario!.(sim_scens[i]["inf_profile_isolation"],
            sim_scens[i]["infection_profiles"],  sim_scens[i]["test_pos_prob"], sim_scens[i]["test_result_days"], 
            sim_scens[i]["symp_day"] .+ 1,  sim_scens[i]["will_isolate"], sim_scens[i]["VL_profiles"], 
            sim_scens[i]["conf_PCR"]; Day7release=Day7release_bool, Day67tests=Day67tests_bool)
    end
    sim_baseline["test_pos_prob"] = fill(zeros(0),sim_baseline["Ntot"])
    sim_baseline["test_result_days"]  = fill(zeros(Int64,0),sim_baseline["Ntot"])
    Conf_PCR_h = fill(zeros(Bool,0),sim_baseline["Ntot"])
    sim_baseline["inf_profile_isolation"] = copy.(sim_baseline["infection_profiles"])
    sim_baseline["isol_days"] = run_testing_scenario!.(sim_baseline["inf_profile_isolation"], 
        sim_baseline["infection_profiles"], sim_baseline["test_pos_prob"],
        sim_baseline["test_result_days"], sim_baseline["symp_day"] .+ 1, 
        sim_baseline["will_isolate"], sim_baseline["VL_profiles"], 
        Conf_PCR_h; Day7release=Day7release_bool, Day67tests=Day67tests_bool)
    sim_scens[Nscens] = copy(sim_baseline)
    for i in 1:(Nscens)
        sim_scens[i]["inf_days"] = zeros(Int64, Ntot)
        for j in 1:Ntot
            sim_scens[i]["inf_days"][j] = sum((sim_scens[i]["inf_profile_isolation"][j] .> 0))
        end
    end
    
    return sim_scens, scen_names
end



function get_no_testing_scenario(Ntot::Int, Pisol::Float64; 
        Day7release_bool::Bool = false, Day67tests_bool::Bool = true)
    sim_baseline = init_VL_and_infectiousness(Ntot, Pisol)
    sim_baseline["test_pos_prob"] = fill(zeros(0),sim_baseline["Ntot"])
    sim_baseline["test_result_days"]  = fill(zeros(Int64,0),sim_baseline["Ntot"])
    Conf_PCR_h = fill(zeros(Bool,0),sim_baseline["Ntot"])
    sim_baseline["inf_profile_isolation"] = copy.(sim_baseline["infection_profiles"])
    sim_baseline["isol_days"] = run_testing_scenario!.(sim_baseline["inf_profile_isolation"], 
        sim_baseline["infection_profiles"], sim_baseline["test_pos_prob"],
        sim_baseline["test_result_days"], sim_baseline["symp_day"] .+ 1, 
        sim_baseline["will_isolate"], sim_baseline["VL_profiles"], 
        Conf_PCR_h; Day7release=Day7release_bool, Day67tests=Day67tests_bool)

    return sim_baseline
end

function run_testing_scenarios_vs_baseline(sim_baseline::Dict, LFD_comply::Float64, 
                                      Conf_PCR::Bool; Day7release_bool::Bool = false,
                            Day67tests_bool::Bool = true, LFD_AllorNone::Bool = false)
    Ntot = sim_baseline["Ntot"]
    Nscens = length(scen_names)
    sim_scens = Array{Dict,1}(undef,Nscens)
    for i in 1:(Nscens-1)
        sim_scens[i] = copy(sim_baseline)
        init_testing_random!(sim_scens[i], Dict("scenario"=>scen_names[i],
                 "comply_prob"=>LFD_comply, "aon_compliance"=>LFD_AllorNone), Conf_PCR)
        
        sim_scens[i]["inf_profile_isolation"] = copy.(sim_scens[i]["infection_profiles"])
        sim_scens[i]["isol_days"] = run_testing_scenario!.(sim_scens[i]["inf_profile_isolation"],
            sim_scens[i]["infection_profiles"],  sim_scens[i]["test_pos_prob"], sim_scens[i]["test_result_days"], 
            sim_scens[i]["symp_day"] .+ 1,  sim_scens[i]["will_isolate"], sim_scens[i]["VL_profiles"], 
            sim_scens[i]["conf_PCR"]; Day7release=Day7release_bool, Day67tests=Day67tests_bool)
    end
    sim_scens[Nscens] = copy(sim_baseline)
    for i in 1:(Nscens)
        sim_scens[i]["inf_days"] = zeros(Int64, Ntot)
        for j in 1:Ntot
            sim_scens[i]["inf_days"][j] = sum((sim_scens[i]["inf_profile_isolation"][j] .> 0))
        end
    end
    
    return sim_scens, scen_names
end