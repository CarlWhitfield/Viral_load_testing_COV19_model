"""
    Alt_inf_models.jl

Written by Carl Whitfield, University of Manchester, 2021

Only to be included by definitions.jl
"""

#=================== Infectivity vs. VL params ===================#

#====> Linear or flat inf models <====#
const inf_dep = 1.3           #Marks 2020 VL vs peak inf dep
const PIsigma = 0.1           #infectivity noise
const PImu = -0.5*PIsigma^2   #relative mean
#Viral load where people stop being infectious
const inf_VL_cutoff = 6.0 
#this scaling means that, on average, the infectivity is 1.0 over the 10-day period
const j0scale = 4.7


#=================================================================#
function get_mean_peak_VL()
    if VL_model == ke_model_no
        return Kmvμ[1]/log(10)
    elseif VL_model == kissler_model_no
        return 7.5
    end
end

VL_ref = get_mean_peak_VL()
function generate_peak_infectivity(peak_VL::Float64)
    if peak_inf_opt == marks_peakinf_opt 
        #assume reference person has average infectivity of 1 (so peak of 2)
        value = j0scale * inf_dep^(peak_VL - VL_ref)
        r = rand(LogNormal(PImu,PIsigma))
        inf = r*value
    elseif peak_inf_opt == unif_peakinf_opt
        inf = 2.0
    end 
        
    return inf
end

function infectivity_linear(PVL::Float64, r::Float64, d::Float64, tp::Float64, 
                          PInf::Float64, T::Int64)
    inf = zeros(T)        
    if PVL > inf_VL_cutoff    
        tinf_start = max(0, tp -  log(10)*(PVL - inf_VL_cutoff)/r)
        tinf_end = tp + log(10)*(PVL - inf_VL_cutoff)/d
        cum_inf = zeros(T+1)
        t_inf = collect(-1:(T-1)) .+ 0.5
        
        #ensure cumulative infectiousness is preserved in interpolation
        cond1 = (t_inf .>= tinf_start) .* (t_inf .<  tp)
        cum_inf[cond1] = 0.5 * PInf .* (t_inf[cond1] .- tinf_start).^2 / 
                               (tp - tinf_start)
        
        cond2 = (t_inf .>= tp) .* (t_inf .<=  tinf_end)
        c_inf_peak = 0.5 * PInf * (tp - tinf_start)
        cum_inf[cond2] = c_inf_peak .+ PInf .* (t_inf[cond2] .- tp) .*
                    (1.0 .- 0.5 .* (t_inf[cond2] .- tp) ./ (tinf_end  - tp))
        cum_inf[t_inf .>  tinf_end] .= 0.5*PInf*(tinf_end - tinf_start)
                
        inf = cum_inf[2:(T+1)] .- cum_inf[1:T]
    end
    
    return inf
end

function infectivity_flat(PVL::Float64, r::Float64, d::Float64, tp::Float64, 
                          PInf::Float64, T::Int64)
    inf = zeros(T)
    
    if PVL > inf_VL_cutoff    
        #if peak viral load high enough to be infectious
        #infectiousness is 1 above VL cutoff
        tinf_start = max(0, tp -  log(10)*(PVL - inf_VL_cutoff)/r)
        tinf_end = tp + log(10)*(PVL - inf_VL_cutoff)/d
        cum_inf = zeros(T+1)
        t_inf = collect(-1:(T-1)) .+ 0.5
           
        #ensure cumulative infectiousness is preserved in interpolation
        cond1 = (t_inf .>= tinf_start) .* (t_inf .<=  tinf_end)
        cum_inf[cond1] = PInf.*(t_inf[cond1] .- tinf_start)
        cum_inf[t_inf .>  tinf_end] .= PInf*(tinf_end - tinf_start)
        
        inf = cum_inf[2:(T+1)] .- cum_inf[1:T]
        #nel = 1:length(t)
        #deleteat!(inf, nel[t .> Int64(ceil(tinf_end))])
    end
    
    return inf
end
    
