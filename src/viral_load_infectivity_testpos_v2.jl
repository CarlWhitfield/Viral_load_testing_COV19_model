"""
    viral_load_infectivity_testpos.jl

Written by Carl Whitfield, University of Manchester, 2021
"""

#============== Options (uncomment and change here) ==============#
#VL_model = ke_model_no         #choice of viral load model (Ke et al is default)
#LFD_model = SC_data_2021         #choice of LFD sensitivity model (social care binned 2021 data is default)
#PCR_sens_max = 0.83            #max PCR sensitivity (Ferretti et al 2021)
#Inf_model = ke_inf_model_no    #infectivity model
#p_asymp = 0.5                  #asymptomatic fraction
#onset_opt = V0_model_opt        #IF using kissler VL model, option for how to generate onset time
#peak_inf_opt = marks_peakinf_opt   #IF using flat or linear inf model, what function is used for peak value
#PCR_TaT_scale = 1.0 
#===============================================================#

include("definitions.jl")

function generate_VL_params(Asymp::Bool)
    if VL_model == kissler_model_no
        return generate_VL_params_Kissler(Asymp)
    elseif VL_model == ke_model_no
        return generate_ke_VL_params()
    end
end

"""
    generate_symp_time(tp::Float64)

Generate and return single time of symptom onset, 
bounded by the times where viral load is at detectable levels

## Arguments: 
`tp` = Time from infection to peak VL (days)

## Returns: 
`Float64` = Symptom onset time (days)

## See also:
`generate_VL_params(Asymp::Bool)`
"""
function generate_symp_time(tp::Float64)
    #Truncate in narrow window, 2 days either side of VL peak
    Strunc = truncated(Gamma(symp_alpha,1.0/symp_beta), 
                       max(0.0, tp - 2.0), tp + 2.0)
    return rand(Strunc)
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


function infectivity(PVL::Float64, r::Float64, d::Float64, tp::Float64, 
                     T::Int64)
    if Inf_model == ke_inf_model_no
        J,h = generate_ke_inf_params()
        inf = infectivity_Ke(PVL, r, d, tp, J, h, T)
    elseif Inf_model == flat_inf_model_no
        PInf = generate_peak_infectivity(PVL)
        #divide by 2 so area under inf curve same as linear case
        inf = infectivity_flat(PVL, r, d, tp, 0.5*PInf, T)  
    elseif Inf_model == linear_inf_model_no
        PInf = generate_peak_infectivity(PVL)
        inf = infectivity_linear(PVL, r, d, tp, PInf, T)
    end
    
    return inf
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
    Asymp = generate_asymptomatic()
    PVL, tp, r, d = generate_VL_params(Asymp)
    
    T = Int64(ceil(tp + log(10)*(PVL - VL_LOD_PCR)/d))
    if T < 1
        print("Strange parameter set generated: PVL = ", PVL, 
              ", d = ", d, ", tp = ", tp, '\n')
        T = 1
    end
    v = zeros(T)
    i = 1:T                              #day index
    t = i .- 1                           #days since infection
    
    v[t .<= tp] = PVL .+  (r/log(10))*(t[t .<= tp] .- tp)
    v[t .> tp] = PVL .- (d/log(10))*(t[t .> tp] .- tp)
    
    ST = generate_symp_time(tp)
    inf = infectivity(PVL, r, d, tp, T)
    SD = Int64(round(ST))

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
               "symp_day"=>zeros(Int64,Ntot),
               "isolation_time"=>zeros(Int64, Ntot),
               "will_isolate"=>zeros(Bool, Ntot),
               "inf_mag"=>zeros(Float64, Ntot),
               "infection_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
               "VL_profiles"=>Array{Array{Float64,1},1}(undef,Ntot),
               "days_infectious" => zeros(Int64,Ntot))
#     for i in 1:Ntot
#         push!(sim["infection_profiles"],zeros(0))
#         push!(sim["VL_profiles"],zeros(0))
#     end
    build_viral_load_distributions!(sim)
    sim["will_isolate"][generate_isolations(Ntot, Pisol)] .= true
    nr = 1:sim["Ntot"]
    sim["non_isolators"] = nr[(sim["will_isolate"] .== false)]
    sim["will_isolate"][sim["asymptomatic"]] .= false #asymptomatics don't self isolate, but are not "non isolators"

    return sim
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
    
function plot_trajectories(A::Array{Array{Float64,1},1})
    Ntraj = length(A) 
    tmax = max(length.(A)...)
    traj = zeros(tmax,Ntraj)
    for n in 1:Ntraj
        traj[1:length(A[n]),n] = A[n]
    end
    npick = collect(1:Ntraj)
    if Ntraj > 100
        npick = rand(npick,100)
    end
    
    plot(1:tmax, traj[:,npick],label=:none, color=:grey, alpha=0.3)
    plot!(1:tmax, median(traj,dims=2), style=:dashdot, color=:black,
          linewidth=3, label="median")
    p1 = plot!(1:tmax, mean(traj,dims=2), style=:dot, color=:red, linewidth=3,
        xlabel="Days since infection", label="mean")
    return p1
end