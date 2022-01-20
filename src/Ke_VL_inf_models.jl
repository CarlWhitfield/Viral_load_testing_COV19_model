"""
    Ke_VL_model.jl

Written by Carl Whitfield, University of Manchester, 2021

Only to be included by definitions.jl
"""

#===================== VL trajectory params  =====================#

#====> Ke et al. 2021 (fitted) <====#
#log means of (1) peak viral load (copies/ml), (2) peak time (days from infection)
#(3) 1/r -- inverse growth rate (days) (4) 1/d -- inverse decay rate (days)
const Kmvμ = [17.546426, 1.392914, -1.197782, -0.671588]   
#log covariance matrix
const KmvΣ = [0.909403 0.043808 -0.010543 0.080628; 
              0.043808 0.032527 0.028866 0.004503; 
             -0.010543 0.028866 0.028640 0.000383; 
              0.080628 0.004503 0.000383 0.026590]
const Log_VLparams_dist = MvNormal(Kmvμ, KmvΣ)  #Distribution used to generate params

#=================================================================#

#=================== Infectivity vs. VL params ===================#

#====> Ke et al. 2021 (fitted) <====#
#log means of parameters J and h (Inf = J VL^h / (VL^h + Km^h))
const Kmvinfμ = [1.782, -0.108]
#log covars
const KmvinfΣ = [0.3387 0.0265; 0.0265 0.1270]
const Jscale = 4.188/exp(Kmvinfμ[1] + KmvinfΣ[1,1]/2)  #normalisation of mean peak inf (a.u.)
const Kinf_Km = 4.0e6
const Log_infparams_dist = MvNormal(Kmvinfμ, KmvinfΣ)  #Distribution used to generate params

# KmvΣ = reshape(KmvΣ,(4,4))
function generate_ke_VL_params()
    Ps = exp.(rand(MvNormal(Kmvμ, KmvΣ)))
    return log10(Ps[1]), Ps[2], 1/Ps[3], 1/Ps[4]
end

function generate_ke_inf_params()
    return exp.(rand(MvNormal(Kmvinfμ, KmvinfΣ)))
end

function infectivity_Ke(PVL::Float64, r::Float64, d::Float64, tp::Float64, J::Float64, h::Float64, T::Int64)
    #if viral load is piecewise linear, and infectivity is hill function, what is cumulative infectivity?
    m_up = r/log(10)
    m_down = d/log(10)

    cum_inf = zeros(T+1)
    t_inf = collect(-1:(T-1)) .+ 0.5
    t_inf[1] = 0
    t_inf[T+1] = T - 1.0
    
    V0 = (PVL - m_up*tp)
    K_m = 4.0e6^h
    a = h*r
    PInf = Jscale*J
    c1 = K_m * 10^(-V0*h)
    cond1 = (t_inf .> 0) .* (t_inf .<=  tp)
    cum_inf[cond1] = PInf .* (log.(exp.(a .* t_inf[cond1]) .+ c1) 
                              .- log(c1 + 1.0)) ./ a
    ci_peak = PInf .* (log.(exp.(a * (tp)) .+ c1) .- log(c1 + 1)) ./ a
    cond2 = (t_inf .>  tp)
    b = -h*d
    c2 = K_m * 10^(-PVL*h)
    cum_inf[cond2] = ci_peak .+ PInf .* (log.(exp.(b .* (t_inf[cond2] .- tp)) .+ c2) 
                              .- log(c2 + 1.0)) ./ b
    
    inf = (cum_inf[2:(T+1)] .- cum_inf[1:T]) ./ (t_inf[2:(T+1)] .- t_inf[1:T])
    
    return inf
end
