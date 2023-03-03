include("./repeat_testing_SC_with_shifts.jl")
using DataFrames
using StatsPlots
using CSV
using Printf
using StatsBase

#key options
VL_model = ke_model_no         #choice of viral load model (Ke et al is default)
LFD_model = porton_down_p3b         #choice of LFD sensitivity model (social care binned 2021 data is default)
PCR_sens_max = 0.95            #max PCR sensitivity (Ferretti et al 2021)
Inf_model = ke_inf_model_no    #infectivity model

function define_parameter_ranges(p::Int)
    NormRange = collect(-0.5:(1.0/(p-1)):0.5)    #normalised parameter range (-0.5:0.5)
    
    #viral load parameter ranges
    VPμrange = log.(10 .^(log10.(exp(DefaultKmvμ[1])) .+ NormRange.*2.4))  #10^(Vp +/- 1.2) 
    TPμrange = log.(exp(DefaultKmvμ[2]) .+ NormRange.*2.0)    #peak time +/- 1 days (fairly well established)
    TGμrange = log.(exp(DefaultKmvμ[3]) .+ NormRange.*0.1)    #1/growth rate +/- 0.05 (changes V0)
    TDμrange = log.(exp(DefaultKmvμ[4]) .+ 0.2 .+ 0.6*NormRange)  #1/decay rate +/- 0.3 (shifted to longer decay times)
    
    #infectiousness parameter ranges
    hμrange = DefaultKmvinfμ[2] .+ NormRange.*2.4      #log(h) +/- 1.2
    Kmμrange = 4.0e6 .* 10 .^(NormRange.*2.4)   #km 4 * 10.^(6 +/- 1.2)
    
    #lfd sensitivity
    VL0Range = DefaultLFD_VL0_pd3b .+ 1 .+ NormRange.*3.0   #cutoff for detection +/- 1.5 log10 copies/ml (shifted for higher cutoff)
    SensRange = DefaultLFD_self_sens .- 0.05 .+ NormRange.*0.3    #max sensitivity +/- 0.15 (shifted for higher cutoff)
    VLkRange = exp.(log(DefaultLFD_VLk_pd3b) .+ NormRange.*1.2)   #steepness of sensitivity slope log(k) +/- 0.6
    
    #other parameters
    SympFracRange = 0.5 .+ NormRange.*0.6         #symptomatic fraction +/- 0.3
    SympMeanRange = 4.84 .+ NormRange.*3.0        #mean onset +/- 1.5 days (subject to constraints)
    SympBetaRange = SympMeanRange ./ 2.6^2
    SympAlphaRange = SympMeanRange .* SympBetaRange
    
    #collect all ranges
    AllLevels = hcat(VPμrange,TPμrange,TGμrange,TDμrange,hμrange,Kmμrange,SensRange,VL0Range,
                 VLkRange,SympFracRange,SympMeanRange)
    varnames = ["Peak VL","Peak VL time","VL growth rate","VL decay rate","inf_h","inf_Km","LFD_sensmax","LFD_sens_VL0",
            "LFD_sens_VLk","Symp_frac","Symp_day"]
    
    return varnames, AllLevels
end

function build_paths(r::Int, p::Int, k::Int, AllLevels::Array{Float64,2})
    #build r paths
    AllPaths = zeros(Int64,k,k+1,r)
    for ir in 1:r
        AllPaths[:,1,ir] = rand(1:p,size(AllLevels)[2])
        order = shuffle(1:k)
        for (n,ik) in enumerate(order)
            AllPaths[:,1+n,ir] = AllPaths[:,n,ir]
            if AllPaths[ik,1,ir] < 1 + p/2
                AllPaths[ik,1+n,ir] = AllPaths[ik,1,ir] + p/2
            else
                AllPaths[ik,1+n,ir] = AllPaths[ik,1,ir] - p/2
            end
        end
    end    
    
    return AllPaths
end

function path_dist(path1::Array{Int64,2},path2::Array{Int64,2})
    #RMSD between paths
    d = 0
    for i in 1:size(path1)[2]
        for j in 1:size(path2)[2]
            d = d + sqrt(sum((path1[:,i] - path2[:,j]).^2))
        end
    end
    return d
end

function path_spread(r::Int, p::Int, AllPaths::Array{Int64,3})
    #measure spread of paths
    d2 = 0
    for irm in 1:r
        for irl in irm:r
            dh = path_dist(AllPaths[:,:,irm],AllPaths[:,:,irl])
            d2 = d2 + (dh/p)^2
        end
    end 
    
    return d2
end

function generate_and_replace_paths!(AllPaths::Array{Int64,3}, m::Int, r::Int, p::Int, k::Int, d2::Float64)
    #generate new paths and replace if they improve the spread of the r paths
    Current_d2 = d2
    for im in 1:m
        NewPath = zeros(Int64,k,k+1)
        NewPath[:,1] = rand(1:p,k)
        order = shuffle(1:k)
        for (n,ik) in enumerate(order)
            NewPath[:,1+n] = NewPath[:,n]
            if NewPath[ik,1] < 1 + p/2
                NewPath[ik,1+n] = NewPath[ik,1] + p/2
            else
                NewPath[ik,1+n] = NewPath[ik,1] - p/2
            end
        end
        New_d2s = Current_d2*ones(r)
        for ir_replace in 1:r
            New_d2s[ir_replace] = Current_d2
            for irm in 1:r
                if ir_replace != irm
                    old_dh = path_dist(AllPaths[:,:,ir_replace],AllPaths[:,:,irm])
                    new_dh = path_dist(NewPath,AllPaths[:,:,irm])
                    New_d2s[ir_replace] = New_d2s[ir_replace] - (old_dh/p)^2 + (new_dh/p)^2
                end
            end
        end
        if max(New_d2s...) > Current_d2
            ir_replace = findmax(New_d2s)[2]
            AllPaths[:,:,ir_replace] = NewPath
            Current_d2 = max(New_d2s...)
        end
    end
    
    return Current_d2
end

function run_sensitivity_sims(r::Int, p::Int, k::Int, Ntot::Int, AllPaths::Array{Int64,3},
                              AllLevels::Array{Float64,2}, varnames::Array{String,1})
    Mean1 = zeros(k+1,r)
    Mean2 = zeros(k+1,r)
    for ir in 1:r
        for ik in 1:(k+1)
            global Kmvμ = [AllLevels[AllPaths[1,ik,ir],1],AllLevels[AllPaths[2,ik,ir],2],
                    AllLevels[AllPaths[3,ik,ir],3],AllLevels[AllPaths[4,ik,ir],4]]
            global Kmvinfμ[2] = AllLevels[AllPaths[5,ik,ir],5]
            global Kinf_Km = AllLevels[AllPaths[6,ik,ir],6]
            global LFD_self_sens = AllLevels[AllPaths[7,ik,ir],7]
            global LFD_VL0_pd3b = AllLevels[AllPaths[8,ik,ir],8]
            global LFD_VLk_pd3b = AllLevels[AllPaths[9,ik,ir],9]
            global Pasymp = 1-AllLevels[AllPaths[10,ik,ir],10]
            SympMeanTime = AllLevels[AllPaths[11,ik,ir],11]
            global symp_beta = SympMeanTime / 2.6^2
            global symp_alpha = SympMeanTime * symp_beta

            sim_baseline = get_no_testing_scenario(Ntot, 1.0)
            scenarios, names = run_testing_scenarios_vs_baseline(sim_baseline, 1.0, false; ScensToRun=[5,6,9])

            tot_IP0, tot_IP1, tot_IP2 = 0, 0, 0
            for i in 1:Ntot
                tot_IP0 = tot_IP0 + sum(scenarios[3]["inf_profile_isolation"][i])
                tot_IP1 = tot_IP1 + sum(scenarios[1]["inf_profile_isolation"][i])
                tot_IP2 = tot_IP2 + sum(scenarios[2]["inf_profile_isolation"][i])
            end
            Mean1[ik,ir] = (tot_IP0 - tot_IP1)/(tot_IP0)
            Mean2[ik,ir] = (tot_IP0 - tot_IP2)/(tot_IP0)
            if isnan(Mean1[ik,ir])
               print("NaN case: VL params: ", Kmvμ, "\nInf params: ", Kmvinfμ, "\nKm: ", Kinf_Km,
                     "\nLFD sens: ", LFD_self_sens, "\nLFD V0: ", LFD_VL0_pd3b, 
                     "\nLFD Vk: ", LFD_VLk_pd3b, "\n Pasymp: ", Pasymp, "\n Symp time: ", 
                     SympMeanTime, '\n') 
            end
        end
    end
    EE1 = zeros(k,r)
    EE2 = zeros(k,r)
    klist = collect(1:k)
    for ir in 1:r
        for ik in 1:k
            kdiff = klist[(AllPaths[:,ik+1,ir] - AllPaths[:,ik,ir]) .!= 0][1]
            ksign = sign(AllPaths[kdiff,ik+1,ir] - AllPaths[kdiff,ik,ir])
            EE1[kdiff,ir] = ksign*(Mean1[ik+1,ir] - Mean1[ik,ir])
            EE2[kdiff,ir] = ksign*(Mean2[ik+1,ir] - Mean2[ik,ir])
        end
    end
    Delta = 1/(p-1)
    EE1 = EE1/(4*Delta)
    EE2 = EE2/(4*Delta)
    mus1 = mean(abs.(EE1),dims=2)
    mus2 = mean(abs.(EE2),dims=2)
    mu1 = mean(EE1,dims=2)
    mu2 = mean(EE2,dims=2)
    sd1 = std(EE1,dims=2)
    sd2 = std(EE2,dims=2)
    
    #print out results to file
    mat1 = hcat(varnames,mus1,mu1,sd1)
    mat2 = hcat(varnames,mus2,mu2,sd2)
    
    return mat1, mat2
end

function run_sensitivity_analysis(repeat::Int, varnames::Array{String,1}, AllLevels::Array{Float64,2},
                                  m::Int, r::Int, p::Int, k::Int, Nsim::Int)
    AllPaths = build_paths(r, p, k, AllLevels)
    print("Paths built\n")
    
    #measure spread of paths
    Current_d2 = path_spread(r, p, AllPaths)
    
    #generate new paths and replace if they improve the spread of the r paths
    Current_d2 = generate_and_replace_paths!(AllPaths, m, r, p, k, Current_d2)
    print("Paths replaced\n")
    
    #run simulations for the r paths and calculate summary stats
    Smat_DeltaIP_2LFDs, Smat_DeltaIP_DailyLFDs = run_sensitivity_sims(r, p, k, Nsim, 
                                                                      AllPaths, AllLevels, varnames)
    print("Sims finished\n")
    
    fname1 = Printf.@sprintf("sens_analysis_2LFDs_repeat%d.csv",repeat)
    fname2 = Printf.@sprintf("sens_analysis_dailyLFD_repeat%d.csv",repeat)
    CSV.write(fname1,Tables.table(Smat_DeltaIP_2LFDs))
    CSV.write(fname2,Tables.table(Smat_DeltaIP_DailyLFDs))
    
    return Smat_DeltaIP_2LFDs, Smat_DeltaIP_DailyLFDs
end

function bootstrap_stderr(data::Array{Float64,1}, B::Int)
    means = zeros(B)
    for b in 1:B
        samples = sample(data,length(data),replace=true)
        means[b] = mean(samples)
    end
    
    return std(means)
end

function run_and_summarise_sensitivity_analyses(Nrepeats::Int64)
    Np = 8           #p in EE analysis (number of levels for each parameter)
    Nrpaths = 40   #r in EE analysis (number of paths to simulate)
    Nmpaths = 2000   #~m in EE analysis (number of new paths to try)
    Ntot = 1000  #number of simulations to run to estimate IP
    varnames, AllLevels = define_parameter_ranges(Np)
    Nparams = size(AllLevels)[2]
    
    repeats = 1:Nrepeats
    outputs = run_sensitivity_analysis.(repeats, Ref(varnames), Ref(AllLevels), Ref(Nmpaths), 
                                        Ref(Nrpaths), Ref(Np), Ref(Nparams), Ref(Ntot))
    k = size(outputs[1][1],1)
    S1mean = zeros(k,3)
    S2mean = zeros(k,3)
    for n in repeats
        S1mean = S1mean .+ outputs[n][1][:,2:4]
        S2mean = S2mean .+ outputs[n][2][:,2:4]
    end
    S1mean = S1mean./Nrepeats
    S2mean = S2mean./Nrepeats
    #bootstrap estimate of stderr
    
    S1stderr = zeros(k,3)
    S2stderr = zeros(k,3)
    for ik in 1:k
        for ij in 1:3
            measured1 = zeros(length(repeats))
            measured2 = zeros(length(repeats))
            for n in repeats
                measured1[n] = outputs[n][1][ik,ij+1]
                measured2[n] = outputs[n][2][ik,ij+1]
            end
            S1stderr[ik,ij] = bootstrap_stderr(measured1, 100)
            S2stderr[ik,ij] = bootstrap_stderr(measured2, 100)
        end
    end
    
    #print output to file
    latex_names = ["Median peak VL \$V_p\$","Median peak VL time \$\\tau_p\$",
    "Median VL inv. growth \$1/r\$", "Median VL inv. decay \$1/d\$",
    "Median inf. sigmoidal slope \$h\$", "Inf. scale param. \$K_m\$", 
    "LFD max. sens. \$\\lambda\$", "LFD sens. cutoff \$\\mu_l\$", 
    "LFD sigmoidal slope \$s_l\$", "Symp. prob. \$P_{\\rm symp}\$",
    "Symp. onset \$\\mu_s\$"]
    
    #do file output here
    fname1 = "latex_output_2LFDs.txt"
    f = open(fname1,"w")
    S1order = reverse(sortperm(vec(S1mean[:,1])))
    for n in S1order
        Printf.@printf(f,"%s & %.3f \\pm %.3f & %.3f \\pm %.3f & %.3f \\pm %.3f\\\\ \\hline \n",     
                    latex_names[n], S1mean[n,1],S1stderr[n,1],S1mean[n,2],S1stderr[n,2],
                    S1mean[n,3],S1stderr[n,3])
    end
    close(f)
    
    fname2 = "latex_output_DailyLFDs.txt"
    f = open(fname2,"w")
    S2order = reverse(sortperm(vec(S2mean[:,1])))
    for n in S2order
        Printf.@printf(f,"%s & %.3f \\pm %.3f & %.3f \\pm %.3f & %.3f \\pm %.3f\\\\ \\hline \n",     
                    latex_names[n], S2mean[n,1],S2stderr[n,1],S2mean[n,2],S2stderr[n,2],
                    S2mean[n,3],S2stderr[n,3])
    end
    close(f)
    
    return S1mean, S1stderr, S2mean, S2stderr
end