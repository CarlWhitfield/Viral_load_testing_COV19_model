"""
    Ke_VL_inf_models.jl

Written by Carl Whitfield, University of Manchester, 2021

Only to be included by definitions.jl
"""

#===================== VL trajectory params  =====================#

#====> Ke et al. 2021 (fitted) <====#
#log means of (1) peak viral load (copies/ml), (2) peak time (days from infection)
#(3) 1/r -- inverse growth rate (days) (4) 1/d -- inverse decay rate (days)
const DefaultKmvμ = ReadOnlyArray([17.546426, 1.392914, -1.197782, -0.671588]) 
Kmvμ = Array(DefaultKmvμ)
#log covariance matrix
const DefaultKmvΣ = Array([0.909403 0.043808 -0.010543 0.080628; 
        0.043808 0.032527 0.028866 0.004503; 
       -0.010543 0.028866 0.028640 0.000383; 
        0.080628 0.004503 0.000383 0.026590])
KmvΣ = deepcopy(DefaultKmvΣ)

#=================================================================#

#=================== Infectivity vs. VL params ===================#

#====> Ke et al. 2021 (fitted) <====#
#log means of parameters J and h (Inf = J VL^h / (VL^h + Km^h))
const DefaultKmvinfμ = [1.782, -0.108]
Kmvinfμ = DefaultKmvinfμ
#log covars
const DefaultKmvinfΣ = [0.3387 0.0265; 0.0265 0.1270]
KmvinfΣ = DefaultKmvinfΣ
const Jhcorr = KmvinfΣ[1,2]/sqrt(KmvinfΣ[1,1]*KmvinfΣ[2,2])
JscaleKe = 5.860/exp(Kmvinfμ[1] + KmvinfΣ[1,1]/2)  #normalisation of mean peak inf (a.u.)
JscaleKiss = 1.090*JscaleKe
const DefaultKinf_Km = 4.0e6
Kinf_Km = DefaultKinf_Km

""" 
    generate_ke_VL_params

Generate parameters for VL model. 

## Returns:
`Float` = peak viral load (log10 copies/ml)
`Float` = peak viral load time (days since infection)
`Float` = growth rate (1/days)
`Float` = decay rate (1/days)
"""
function generate_ke_VL_params()
    Ps = exp.(rand(MvNormal(Kmvμ, KmvΣ)))
    return log10(Ps[1]), Ps[2], 1/Ps[3], 1/Ps[4]
end

""" 
    generate_ke_inf_params_correlated

Generate parameters for infectiousness model assuming 1:1 correlation between
peak viral load and peak infectivity. 

## Arguments
`PVL::Float64` = Peak viral load (log10 copies/ml)
`PVLμ::Float64` = Mean of peak viral load distribution
`PVLσ::Float64` = Std of peak viral load distribution

## Returns:
`Float` = J_p (infectiousness scale)
`Float` = h (Hill parameter)
"""
function generate_ke_inf_params_correlated(PVL::Float64, PVLμ::Float64, PVLσ::Float64)
    lJ = Kmvinfμ[1] + (PVL - PVLμ)*sqrt(KmvinfΣ[1,1])/PVLσ
    RandNorm = rand(Normal(0.0,1.0))
    lh = Kmvinfμ[2] + Jhcorr*(lJ - Kmvinfμ[1])*sqrt(KmvinfΣ[2,2]/KmvinfΣ[1,1]) + 
                          sqrt(1-Jhcorr^2)*RandNorm*sqrt(KmvinfΣ[2,2])
    return exp.([lJ,lh])
end

""" 
    generate_ke_inf_params_uncorrelated

Generate infectivity based on derived distribution

## Returns:
`Float` = J_p (infectiousness scale)
`Float` = h (Hill parameter)
"""
function generate_ke_inf_params_uncorrelated()
    return exp.(rand(MvNormal(Kmvinfμ, KmvinfΣ)))
end

""" 
    infectivity_Ke

Generate parameters for infectiousness model assuming 1:1 correlation between
peak viral load and peak infectivity. 

## Arguments
`PVL::Float64` = Peak viral load (log10 copies/ml)
`r::Float64` = growth rate (1/days)
`d::Float64` = decay rate (1/days)
`tp::Float64` = Time to peak viral load (days)
`J::Float64` = infectiousness scale
`h::Float64` = Hill parameter
`T::Int64` = Max number of days to fill
`scale_opt::Int` [optional keyword arg.] = 
     (0: scale for Ke VL model, 1 scale for Kissler VL model)

## Returns:
`Float` = J_p (infectiousness scale)
`Float` = h (Hill parameter)
"""
#if viral load is piecewise exponential, and infectivity is hill function
#so we can integrate a single day timestep analytically
function infectivity_Ke(PVL::Float64, r::Float64, d::Float64, tp::Float64, J::Float64, h::Float64, T::Int64; scale_opt::Int = 0)
    #
    m_up = r/log(10)   #VL linear slope up (log10 copies/ml scale)
    m_down = d/log(10)   #VL linear slope down (log10 copies/ml scale)

    cum_inf = zeros(T+1)   #init cumulative infectiousness
    t_inf = collect(-1:(T-1)) .+ 0.5  #days since infection + 1/2
    t_inf[1] = 0
    t_inf[T+1] = T - 1.0  #days since infection - 1/2
    
    V0 = (PVL - m_up*tp)     #V at t = 0
    K_m = Kinf_Km^h         #scale parameter (Ke et al 2021)
    a = h*r
    if scale_opt == 0
        PInf = JscaleKe*J
    else   
        PInf = JscaleKiss*J
    end
    #calculate integral
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
    
    #daily infectiousness is gradient of integral over the 24h period
    inf = (cum_inf[2:(T+1)] .- cum_inf[1:T]) ./ (t_inf[2:(T+1)] .- t_inf[1:T])
    
    return inf
end
