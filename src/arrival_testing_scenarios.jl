include("viral_load_infectivity_testpos.jl")

const arrival_scen_names = ["Day 2 PCR","Daily LFD","Both"]

VL_model = kissler_model_no 
Inf_model = flat_inf_model_no
peak_inf_opt = unif_peakinf_opt

function arrival_scenario_1_setup(start_day::Int, max_PCR_day::Int, pre_LFD::Bool = false,
                                  pre_LFD_time::Int = 2)   #1 PCR in first 2 days
    TestDays = [rand(start_day:(start_day+max_PCR_day))]
    PCR_delays = draw_PCR_delays(1)
    ResultDays = TestDays .+ PCR_delays
    TestTypes = [0]
    PreIsolationDays = collect(start_day:(ResultDays[1]-1))
    if pre_LFD
        if start_day > pre_LFD_time
            push!(TestDays, start_day-pre_LFD_time)
            push!(ResultDays, start_day-pre_LFD_time)
            push!(TestTypes,1)
        end
    end
    itd = sortperm(ResultDays)
    
    return TestDays[itd], TestTypes[itd], ResultDays[itd], PreIsolationDays
end

function arrival_scenario_2_setup(start_day::Int, NLFDs::Int, pre_LFD::Bool = false,
                                  pre_LFD_time::Int = 2)   #Daily LFDs
    TestDays = collect(start_day:(start_day+NLFDs-1))
    ResultDays = copy(TestDays)
    TestTypes = ones(Int64,NLFDs)
    PreIsolationDays = []
    if pre_LFD
        if start_day > pre_LFD_time
            push!(TestDays, start_day-pre_LFD_time)
            push!(ResultDays, start_day-pre_LFD_time)
            push!(TestTypes,1)
        end
    end
    itd = sortperm(ResultDays)
    
    return TestDays[itd], TestTypes[itd], ResultDays[itd],
           PreIsolationDays
end

function arrival_scenario_3_setup(start_day::Int, NLFDs::Int, max_PCR_day::Int, 
                                 pre_LFD::Bool = false, pre_LFD_time::Int = 2,
                                 LFD_after_PCR::Bool = false)  #Daily LFDs and 1 PCR in first 2 days
    PCR_test_days = [rand(start_day:(start_day+max_PCR_day))]
    PCR_delays = draw_PCR_delays(length(PCR_test_days))
    PCR_result_days = PCR_test_days .+ PCR_delays
    PreIsolationDays = collect(start_day:(PCR_result_days[1]-1))
    
    LFD_test_days = collect(start_day:(start_day+NLFDs-1))  
    if LFD_after_PCR
        LFD_test_days = LFD_test_days .+ PCR_result_days[1]
    end
    if pre_LFD
        if start_day > pre_LFD_time
            push!(LFD_test_days, start_day-pre_LFD_time)
        end
    end
    TestDays = vcat(LFD_test_days,PCR_test_days)
    ResultDays = vcat(LFD_test_days,PCR_result_days)
    TestTypes = vcat(ones(Int64,length(LFD_test_days)), zeros(Int64,length(PCR_test_days)))
    itd = sortperm(ResultDays)
    
    return TestDays[itd], TestTypes[itd], ResultDays[itd],
           PreIsolationDays 
end

function init_testing_random_start!(sim::Dict, testing_params::Dict; 
                                    Ndaysel::Int = 14)
    #remove people who are either not infectious or would be isolating
    #ToKeep = ones(Bool,sim["Ntot"])
    sim["will_isolate_with_test"] = ones(Bool,sim["Ntot"])
    sim["testing_paused"] = zeros(Bool,sim["Ntot"])
    sim["resume_testing"] = -ones(Int64,sim["Ntot"])
    sim["conf_PCR"] = Array{Array{Bool,1},1}(undef,sim["Ntot"])
    sim["test_days"] = Array{Array{Int64,1},1}(undef,sim["Ntot"])
    sim["test_types"] = Array{Array{Int64,1},1}(undef,sim["Ntot"])
    sim["test_result_days"] = Array{Array{Int64,1},1}(undef,sim["Ntot"])
    sim["preisolation"] = Array{Array{Int64,1},1}(undef,sim["Ntot"])
    TestOutput = []
    
    if testing_params["scenario"] == arrival_scen_names[1]
        if testing_params["pre_LFD"]
            TestOutput = arrival_scenario_1_setup.(sim["start_days"], 
                  testing_params["PCR_max_day"]*ones(Int64,sim["Ntot"]), true, 
                  testing_params["pre_LFD_day"])
        else
            TestOutput = arrival_scenario_1_setup.(sim["start_days"], 
                  testing_params["PCR_max_day"]*ones(Int64,sim["Ntot"]))
        end
    elseif testing_params["scenario"] == arrival_scen_names[2]
        if testing_params["pre_LFD"]
            TestOutput = arrival_scenario_2_setup.(sim["start_days"], 
                testing_params["LFD_tests"]*ones(Int64,sim["Ntot"]), true, 
                testing_params["pre_LFD_day"])
        else
            TestOutput = arrival_scenario_2_setup.(sim["start_days"], 
                testing_params["LFD_tests"]*ones(Int64,sim["Ntot"]))
        end
    elseif testing_params["scenario"] == arrival_scen_names[3]
        if testing_params["pre_LFD"]
            TestOutput = arrival_scenario_3_setup.(sim["start_days"], 
                testing_params["LFD_tests"]*ones(Int64,sim["Ntot"]),
                testing_params["PCR_max_day"]*ones(Int64,sim["Ntot"]), Ref(true), 
                Ref(testing_params["pre_LFD_day"]), Ref(testing_params["LFD_after_PCR"]))
        else
            TestOutput = arrival_scenario_3_setup.(sim["start_days"], 
                testing_params["LFD_tests"]*ones(Int64,sim["Ntot"]),
                testing_params["PCR_max_day"]*ones(Int64,sim["Ntot"]), Ref(false), 
                Ref(0), Ref(testing_params["LFD_after_PCR"]))
        end
    end
    
    for i in 1:sim["Ntot"]
        sim["test_days"][i] = TestOutput[i][1]
        sim["test_types"][i] = TestOutput[i][2]
        sim["test_result_days"][i] = TestOutput[i][3]
        sim["conf_PCR"][i] = zeros(Bool,length(sim["test_days"][i]))
#         if Conf_PCR
#             sim["conf_PCR"][i][sim["test_types"][i] .== 1] .= true
#         end
        sim["preisolation"][i] = TestOutput[i][4]
    end
    
    sim["test_pos_prob"] = get_test_probs.(sim["VL_profiles"], sim["test_days"], 
        sim["test_types"]; Do_LFD_test_prob = testing_params["comply_prob"]) 
    #assuming random compliance
end


function run_testing_flight_arrivals(Ntot::Int64, PIsol::Float64, 
             testing_params::Array{Dict,1}, arrival_GR::Float64, Ndaysel::Int,
             infected_in_flight::Bool=false)
    
    sim_baseline = init_VL_and_infectiousness(Ntot, PIsol)
    VL_len = length.(sim_baseline["VL_profiles"])
    #define start day (a.k.a arrival day)
    #exponential distribution (growth rate arrival_GR) -- 0 if uniform
    if infected_in_flight
        sim_baseline["start_days"] = ones(Int64,sim_baseline["Ntot"])
    else
        sim_baseline["start_days"] = sample(2:(Ndaysel+1),Weights(exp.(-arrival_GR.*collect(0:13))),
                                        sim_baseline["Ntot"])
    end
    sim_baseline["days_rel_arrival"] = Array{Array{Int64,1},1}(undef,sim_baseline["Ntot"])
    for i in 1:sim_baseline["Ntot"]
        sim_baseline["days_rel_arrival"][i] = collect(1:VL_len[i]) .- sim_baseline["start_days"][i]
    end
    sims = Array{Dict,1}(undef,length(testing_params))
    for (n,tp) in enumerate(testing_params)
        sim_scen = copy(sim_baseline)
        init_testing_random_start!(sim_scen, tp)
        #Need to sort out symp_day
        
        isolation_days = Array{Array{Int64,1},1}(undef,Ntot)
        for i in 1:Ntot
            isolation_days[i] = zeros(0)
        end
        sim_scen["isol_days"] = run_testing_scenario!.(isolation_days, sim_scen["infection_profiles"], 
            sim_scen["test_pos_prob"], sim_scen["test_result_days"], sim_scen["symp_day"] .+ 1, 
            sim_scen["will_isolate"], sim_scen["VL_profiles"],  sim_scen["conf_PCR"],
            sim_scen["preisolation"])
        
        sim_scen["inf_profile_isolation"] = deepcopy(sim_scen["infection_profiles"])
        for i in 1:Ntot
            valid_isol_days = isolation_days[i][isolation_days[i] .<= length(sim_scen["inf_profile_isolation"][i])]
            sim_scen["inf_profile_isolation"][i][valid_isol_days] .= 0
        end
    
        sims[n] = sim_scen
    end
        
    return sim_baseline, sims
end