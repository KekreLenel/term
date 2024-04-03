# find stochastic SS
function calcSS(prms::params{R,I}, pitp, sitp, aEitp, sigmaitp) where {R <: Real, I <: Integer}


    @unpack  W_center, r_bar, β_bar, Δτ, Δτ_y, τ_y, mu, αvec, n_dims,
             θ0vec, θvec, τ, y_idxs, n_τ, T_ss = prms

    # starting point
    ss = [r_bar, β_bar, W_center]

    # shocks 
    rshocks = zeros(100-1)
    βshocks = zeros(100-1)
    Wshocks = zeros(100-1)

    dif = 1.0 
    while dif > sqrt(eps())
    
        # procede by addding 100 steps of simulation with zero shocks until convergence
        r_series, β_series , W_series = calcSeries(prms, ss, rshocks, βshocks, Wshocks, pitp, sitp, aEitp, sigmaitp, false)

        rss = r_series[end]
        βss = β_series[end]
        Wss = W_series[end]
        ss_new = [rss, βss, Wss]
        dif = sum(abs.(ss_new - ss))
        ss = deepcopy(ss_new)

    end
    
    # one more calculation to extract all relevant series
    r_series, β_series , W_series, y_series, aE_series = calcSeries(prms, ss, rshocks, βshocks, Wshocks, pitp, sitp, aEitp, sigmaitp, false)
    
    # store stochastic steady state
    yss = y_series[y_idxs,end]
    Ess = aE_series[:,1,end]
    fss = (τ_y[2:end].*yss[2:end] - τ_y[1:end-1].*yss[1:end-1])/Δτ_y

    # find wealth shock necessary to stay at W_center with no shocks 
    init = [r_bar, β_bar, W_center]
    dif = 1.0 
    while dif > 1e-12

        r_series, β_series , W_series = calcSeries(prms, init, rshocks, βshocks, Wshocks, pitp, sitp, aEitp, sigmaitp, false)

        dif     = (abs.(W_series[end] - W_center))
        Wshocks = Wshocks .- (W_series[2] - W_center)

    end

    tgt_Wshock = Wshocks[1]

    return Stochsteady(ss, yss, fss, Ess, tgt_Wshock)
    
end


function simSeries(prms::params{R,I}, pitp, sitp, aEitp, sigmaitp, Wss)  where {R <: Real, I <: Integer}

    @unpack  σ_r_unc, T_sim, T_burn, αvec, θvec, r_bar, β_bar, W_bar, W_center, y_idxs, n_idxs, mu,
              τ_y, Δτ, n_τ, τ, Δt_y, Len_sim, N_sim, NN_sim, θ0vec, f_τ, κ_r= prms
    
    # draw simulation shocks
    rng    = MersenneTwister(1234)
    if NN_sim == 1
        sim_draws = [randn(rng, 2, T_sim-1)]
    else
        sim_draws = [randn(rng, 2, T_sim-1) for _ in 1:div(NN_sim,2)]
        sim_draws = vcat(sim_draws, -sim_draws)
    end
    Wshocks = zeros(T_sim-1)

    # initialize series
    moments_vec = [zeros(123) for _ in 1:NN_sim]
    y_avrg_vec  = [zeros(30) for _ in 1:NN_sim]
    fwd_avrg_vec = [zeros(29) for _ in 1:NN_sim]
    β_FB_vec    = [zeros(29) for _ in 1:NN_sim]
    β_CS_vec    = [zeros(29) for _ in 1:NN_sim]
    r_series_all = zeros(T_sim - T_burn, NN_sim)
    β_series_all = zeros(T_sim - T_burn, NN_sim)
    W_series_all = zeros(T_sim - T_burn, NN_sim)
    CP_vec = [zeros(29,30) for _ in 1:NN_sim]
    CPE_vec = [zeros(29,30) for _ in 1:NN_sim]
    
    # loop over simulations (would be more efficient to use @distributed here)
    Threads.@threads for nnn = 1:NN_sim

        rshocks = sim_draws[nnn][1,:]
        βshocks = sim_draws[nnn][2,:]

        # draw long simulation
        init = [r_bar, β_bar, W_center]
        r_series, β_series , W_series, y_series, aE_series, sigma_series = calcSeries(prms, init, rshocks, βshocks, Wshocks, pitp, sitp, aEitp, sigmaitp, true)  

        yy_series = y_series[y_idxs,:] # store annual prices
        
        # get forward and return series 
        fwd_series, ret_series, mret_series = calcRetFwd(prms, y_series, yy_series)

        # burn beginning
        y_series     = y_series[:, T_burn + 1:end]
        yy_series    = yy_series[:, T_burn + 1:end]
        ret_series   = ret_series[:, T_burn + 1:end]
        fwd_series   = fwd_series[:, T_burn + 1:end]
        mret_series  = mret_series[:, T_burn + 1:end]
        r_series     = r_series[T_burn + 1:end]
        β_series     = β_series[T_burn + 1:end]
        W_series     = W_series[T_burn + 1:end]
        rshocks      = rshocks[T_burn + 1:end]
        βshocks      = βshocks[T_burn + 1:end]
        sigma_series = sigma_series[:, T_burn + 1:end]
        
        # store state series for histograms
        r_series_all[:,nnn] .= r_series
        β_series_all[:,nnn] .= β_series
        W_series_all[:,nnn] .= W_series

        # calc regressions moments
        β_FB_vec[nnn][:] .= FBReg(prms, yy_series, fwd_series, ret_series)
        β_CS_vec[nnn][:] .= CSReg(prms, y_series)
        
        # calculate CP decomposition
        CP_tmp   = CPDecomp(prms, y_series, fwd_series, ret_series)
        CP_vec[nnn][:,:]  .= mean(CP_tmp, dims=3)
        CP_tmp   = CPEDecomp(prms, fwd_series, aE_series)
        CPE_vec[nnn][:,:] .= mean(CP_tmp, dims=3)

        # preallocate moments
        moments = zeros(123)

        y_avrg_vec[nnn][:]  .= mean(yy_series, dims=2)[:,1]
        fwd_avrg_vec[nnn][:] .= mean(fwd_series, dims=2)[:,1]

        # yield moments
        
        tmp_mat = reshape(yy_series,n_idxs,Len_sim,N_sim)
        @views tmp1_mat = tmp_mat[1,:,:]

        tmp_idx = 1
        @views moments[1]  = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                                    # 1y yield mean
        @views moments[2]  = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                                    # 1y yield vol
        @views moments[3]  = 100*mean(std( tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[4]  = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])                  # 1y yield yoy ac 
        @views moments[5]  = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])           # avrg corr w 1y

        tmp_idx = 2
        @views moments[6]  = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[7]  = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[8]  = 100*mean(std( tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[9]  = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[10] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_idx = 3
        @views moments[11] = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[12] = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[13] = 100*mean(std( tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[14] = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[15] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_idx = 5
        @views moments[16] = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[17] = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[18] = 100*mean(std( tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[19] = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[20] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_idx = 10
        @views moments[21] = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[22] = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[23] = 100*mean(std( tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[24] = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[25] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_idx = 20
        @views moments[26] = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[27] = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[28] = 100*mean(std( tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[29] = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[30] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        @views Dtmp_mat    = tmp_mat[10,:,:] - tmp_mat[1,:,:]
        @views moments[31] = 100*mean(mean(Dtmp_mat, dims = 1))                               # 1y yield mean
        @views moments[32] = 100*mean(std( Dtmp_mat, dims = 1))                               # 1y yield vol
        @views moments[33] = 100*mean(std( Dtmp_mat[13:end,:] .- Dtmp_mat[1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[34] = mean([autocor(Dtmp_mat[:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[35] = mean([cor(Dtmp_mat[:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        @views Dtmp_mat    = tmp_mat[10,:,:] - tmp_mat[2,:,:]
        @views moments[36] = 100*mean(mean(Dtmp_mat, dims = 1))                               # 1y yield mean
        @views moments[37] = 100*mean(std( Dtmp_mat, dims = 1))                               # 1y yield vol
        @views moments[38] = 100*mean(std( Dtmp_mat[13:end,:] .- Dtmp_mat[1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[39] = mean([autocor(Dtmp_mat[:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[40] = mean([cor(Dtmp_mat[:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        @views Dtmp_mat    = tmp_mat[20,:,:] - tmp_mat[1,:,:]
        @views moments[41] = 100*mean(mean(Dtmp_mat, dims = 1))                               # 1y yield mean
        @views moments[42] = 100*mean(std( Dtmp_mat, dims = 1))                               # 1y yield vol
        @views moments[43] = 100*mean(std( Dtmp_mat[13:end,:] .- Dtmp_mat[1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[44] = mean([autocor(Dtmp_mat[:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[45] = mean([cor(Dtmp_mat[:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        @views Dtmp_mat    = tmp_mat[20,:,:] - tmp_mat[2,:,:]
        @views moments[46] = 100*mean(mean(Dtmp_mat, dims = 1))                               # 1y yield mean
        @views moments[47] = 100*mean(std( Dtmp_mat, dims = 1))                               # 1y yield vol
        @views moments[48] = 100*mean(std( Dtmp_mat[13:end,:] .- Dtmp_mat[1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[49] = mean([autocor(Dtmp_mat[:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[50] = mean([cor(Dtmp_mat[:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_mat = reshape(fwd_series,n_idxs-1,Len_sim,N_sim)

        tmp_idx = 1
        @views moments[51] = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[52] = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[53] = 100*mean(std( tmp_mat[tmp_idx,13:end,:] - tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[54] = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[55] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_idx = 4
        @views moments[56] = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[57] = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[58] = 100*mean(std( tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[59] = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[60] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_idx = 9
        @views moments[61] = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[62] = 100*mean(std( tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[63] = 100*mean(std( tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[64] = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[65] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_idx = 19
        @views moments[66] = 100*mean(mean(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield mean
        @views moments[67] = 100*mean(std(tmp_mat[tmp_idx,:,:], dims = 1))                               # 1y yield vol
        @views moments[68] = 100*mean(std(tmp_mat[tmp_idx,13:end,:] .- tmp_mat[tmp_idx,1:end-12,:], dims = 1)) # d1y yield vol 
        @views moments[69] = mean([autocor(tmp_mat[tmp_idx,:,iii],[1])[1] for iii = 1:N_sim])               # 1y yield yoy ac 
        @views moments[70] = mean([cor(tmp_mat[tmp_idx,:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])  # avrg corr w 1y

        tmp_mat     = reshape(y_series,n_τ,Len_sim,N_sim)
        moments[71] = 100*mean(std(tmp_mat[2:y_idxs[20],:,:],dims=2))  # avrg yield vol 
        moments[72] = 100*mean(std(tmp_mat[2:y_idxs[20],13:end,:].-tmp_mat[2:y_idxs[20],1:end-12,:],dims=2))      # avrg dy vol
        moments[73] = mean([cor(tmp_mat[x,13:end,zzz]-tmp_mat[x,1:end-12,zzz],                 # avrg corr w 1y
                                tmp_mat[y_idxs[1],13:end,zzz]-tmp_mat[y_idxs[1],1:end-12,zzz]) for x = y_idxs[1]:y_idxs[20] for zzz = 1:N_sim])
        

        # trade vol moments
        @views tradevol = Δτ*sum(abs.(-(-αvec.*y_series[:,2:end].*τ   .+ adj_fun.(W_series[2:end]',Ref(prms)).*θvec.*exp.(β_series[2:end]')) 
                                +(-αvec.*y_series[:,1:end-1].*τ .+ adj_fun.(W_series[1:end-1]',Ref(prms)).*θvec.*exp.(β_series[1:end-1]'))), dims=2)
        @views moments[74] = sum(tradevol[1:y_idxs[2]])/sum(tradevol)            # relative 0-2 trade vol 
        @views moments[75] = sum(tradevol[y_idxs[10]+1:y_idxs[end]])/sum(tradevol) # relative 10+ trade vol 
        
        # wealth process moments 
        moments[76] = mean(W_series)                                          #  W mean
        moments[77] = std(log.(W_series))                                     #  W vol
        moments[78] = maximum(log.(W_series))-log.(W_center)                     #  
        moments[79] = minimum(log.(W_series))-log.(W_center)                    #  
        @views moments[80] = cor(W_series[13:end], W_series[1:end-12])               #  dW vol
        moments[81] = Wss

        # position moments
        xmat = -(αvec.*y_series.*τ) .+ θ0vec .+ adj_fun.(W_series',Ref(prms)).*θvec.*β_series'
        @views moments[82] = mean(sum(xmat[1:y_idxs[2],:], dims=1)./sum(xmat, dims=1))
        @views moments[83] = mean(sum(xmat[y_idxs[10]:end,:], dims=1)./sum(xmat, dims=1))
        @views moments[84] = mean(sum(xmat[1:y_idxs[2],:], dims=1))
        @views moments[85] = mean(sum(xmat[y_idxs[10]:end,:], dims=1))

        # duration moments 
        weighted_assets = sum(Δτ.*xmat.*τ,dims=1)
        total_assets    = sum(Δτ.*xmat,dims=1)
        dura            = weighted_assets./total_assets
        dure            = weighted_assets./W_series'

        moments[86]     = mean(dura) 
        moments[87]     = median(dura) 
        moments[88]     = mean(dure) 
        moments[89]     = median(dure)

        # IG
        Amat            = zeros(size(xmat,1), size(xmat,2))
        Amat[xmat.>=0]  = xmat[xmat.>=0]
        Avec            = vec(sum(Amat.*Δτ, dims=1)' .+ (W_series .- sum(xmat.*Δτ, dims=1)').*(W_series .- sum(xmat.*Δτ, dims=1)' .> 0))
        xmat1           = xmat[y_idxs[1]+1:end,:]
        IG              = vec((W_series .- sum(xmat1.*Δτ, dims=1)')./Avec)

        # term premia
        yr5f5           = (10*y_series[y_idxs[10],:] .- 5*y_series[y_idxs[5],:])./5
        tp5f5           = yr5f5 .- ((r_series .- r_bar).*((exp.(-κ_r.*5) .- exp.(-κ_r.*10))./κ_r) .+ r_bar.*5)./5
        tp10            = y_series[y_idxs[10],:] - ((r_series .- r_bar).*((1 .- exp.(-κ_r.*10))./κ_r) .+ r_bar.*10)./10

        # volatility moments

        @views moments[90] = 100*std(log.(W_series[2:end]).-log.(W_series[1:end-1]))
        @views moments[91] = 100*std(log.(W_series[13:end]).-log.(W_series[1:end-12]))
        moments[92]        = 100*std(tp5f5)
        moments[93]        = 100*std(tp10)
        moments[94]        = std(log.(dura[dura.>0]))
        moments[95]        = std(log.(dure[dure.>0]))
        moments[96]        = std(IG)

        tmp_mat            = reshape(log.(W_series),Len_sim,N_sim)
        @views moments[97] = 100*mean(std( tmp_mat[2:end,:] .- tmp_mat[1:end-1,:], dims = 1))
        @views moments[98] = 100*mean(std( tmp_mat[13:end,:] .- tmp_mat[1:end-12,:], dims = 1))

        tmp_mat            = reshape(tp5f5,Len_sim,N_sim)
        moments[99]        = 100*mean(std( tmp_mat, dims = 1))

        tmp_mat            = reshape(tp10,Len_sim,N_sim)
        moments[100]       = 100*mean(std( tmp_mat, dims = 1))

        tmp_mat            = log.(dura[dura.>0])
        max_length         = Int(floor(size(tmp_mat,1)/Len_sim))
        tmp_mat            = reshape(tmp_mat[1:Int(max_length*Len_sim)], Len_sim, max_length)
        moments[101]       = mean(std( tmp_mat, dims = 1))

        tmp_mat            = log.(dure[dure.>0])
        max_length         = Int(floor(size(tmp_mat,1)/Len_sim))
        tmp_mat            = reshape(tmp_mat[1:Int(max_length*Len_sim)], Len_sim, max_length)
        moments[102]       = mean(std( tmp_mat, dims = 1))

        tmp_mat            = reshape(IG,Len_sim,N_sim)
        moments[103]       = mean(std( tmp_mat, dims = 1))

        # conditional volatility moments

        @views moments[104] = 100*std(sigma_series[10,:])
        @views moments[105] = 100*std(sigma_series[20,:])

        @views tmp_mat      = reshape(sigma_series[10,:],Len_sim,N_sim)
        moments[106]        = 100*mean(std( tmp_mat, dims = 1))

        @views tmp_mat      = reshape(sigma_series[20,:],Len_sim,N_sim)
        moments[107]        = 100*mean(std( tmp_mat, dims = 1))

        # term premia moments

        tmp_mat             = reshape(tp5f5,Len_sim,N_sim)
        @views moments[108] = 100*mean(mean(tmp_mat, dims = 1))
        @views moments[109] = 100*mean(std( tmp_mat[13:end,:] .- tmp_mat[1:end-12,:], dims = 1))
        @views moments[110] = mean([autocor(tmp_mat[:,iii],[1])[1] for iii = 1:N_sim])
        @views moments[111] = mean([cor(tmp_mat[:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])

        tmp_mat             = reshape(tp10,Len_sim,N_sim)
        moments[112]        = 100*mean(mean(tmp_mat, dims = 1))
        @views moments[113] = 100*mean(std( tmp_mat[13:end,:] .- tmp_mat[1:end-12,:], dims = 1))
        @views moments[114] = mean([autocor(tmp_mat[:,iii],[1])[1] for iii = 1:N_sim])
        @views moments[115] = mean([cor(tmp_mat[:,iii], tmp1_mat[:,iii]) for iii = 1:N_sim])

        # Yield moments throughout the sample - used for analytical solution
        @views moments[116] = 100*mean(y_series[y_idxs[1],:])
        @views moments[117] = 100*mean(y_series[y_idxs[10],:] .- y_series[y_idxs[1],:])
        @views moments[118] = 100*mean(y_series[y_idxs[20],:] .- y_series[y_idxs[1],:])
        @views moments[119] = 100*std(y_series[y_idxs[5],:])
        @views moments[120] = 100*std(y_series[y_idxs[10],:])
        @views moments[121] = 100*std(y_series[y_idxs[20],:])
        @views moments[122] = 100*std(y_series[y_idxs[20],(f_τ+1):end]-y_series[y_idxs[20],1:(end-f_τ)])

        moments_vec[nnn] = moments
    end 
    
    avrg_moments = mean(moments_vec,dims=1)[1]
    avrg_β_FB    = mean(β_FB_vec,dims=1)[1]
    avrg_β_CS    = mean(β_CS_vec,dims=1)[1]
    avrg_CP      = mean(CP_vec,dims=1)[1]
    avrg_CPE     = mean(CPE_vec,dims=1)[1]
    avrg_y_avrg = mean(y_avrg_vec, dims=1)[1]
    avrg_fwd_avrg = mean(fwd_avrg_vec, dims=1)[1]
    r_series_all = vec(r_series_all)
    β_series_all = vec(β_series_all)
    W_series_all = vec(W_series_all)

    return Simulation(avrg_moments, avrg_β_FB, avrg_β_CS, avrg_y_avrg, avrg_fwd_avrg, r_series_all, β_series_all, W_series_all, avrg_CP, avrg_CPE)

end

function simSeriesQ(prms::params{R,I}, sitp, calib)   where {R <: Real, I <: Integer}

    @unpack  σ_r_unc, T_sim, T_burn, αvec, θvec, r_bar, β_bar, W_bar, W_center, y_idxs, n_idxs, mu,
              τ_y, Δτ, n_τ, τ, Δt_y, Len_sim, N_sim, θ0vec, f_τ, κ_r, κ_β, σ_r, σ_β, τ_max, rbnds,
              NN_sim, ξ, rbnds_x, βbnds, βbnds_x, Wbnds, Wbnds_x, endo_wealth, n_dims = prms
    
    # draw simulation shocks
    rng = MersenneTwister(1234)
    if NN_sim == 1
        sim_draws = [randn(rng, 2, T_sim-1)]
    else
        sim_draws = [randn(rng, 2, T_sim-1) for _ in 1:div(NN_sim,2)]
        sim_draws = vcat(sim_draws, -sim_draws)
    end

    r_series_all = zeros(T_sim - T_burn, NN_sim)
    β_series_all = zeros(T_sim - T_burn, NN_sim)
    W_series_all = zeros(T_sim - T_burn, NN_sim)

    Threads.@threads for nnn = 1:NN_sim

        rshocks = sim_draws[nnn][1,:]
        βshocks = sim_draws[nnn][2,:]

        # preallocation
        r_series  = zeros(T_sim)             # path of interest rates
        β_series  = zeros(T_sim)             # path of demand
        W_series  = zeros(T_sim)             # path of wealth

        # starting point
        r_series[1]  = r_bar
        β_series[1]  = β_bar
        W_series[1]  = W_center

        for ttt = 2:T_sim
            
            if n_dims == 3 
                drifts = sitp(r_series[ttt-1],β_series[ttt-1],log(W_series[ttt-1]))
            else 
                drifts = sitp(r_series[ttt-1],β_series[ttt-1])
            end

            drift_r = drifts.myvec[1] 
            drift_β = drifts.myvec[2] 
            drift_W = drifts.myvec[3] 
            Ereturn = drifts.myvec[4] 
            η_r     = drifts.myvec[5] 
            η_β     = drifts.myvec[6] 

            r_series[ttt] = r_series[ttt-1] + (κ_r*(r_bar-r_series[ttt-1])-drift_r)*Δτ + sqrt(Δτ)*σ_r*rshocks[ttt-1]
            r_series[ttt] = maximum([minimum([r_series[ttt], rbnds_x[2]]), rbnds_x[1]])
            β_series[ttt] = β_series[ttt-1] + (κ_β*(β_bar-β_series[ttt-1])-drift_β)*Δτ + sqrt(Δτ)*σ_β*βshocks[ttt-1]
            β_series[ttt] = maximum([minimum([β_series[ttt], βbnds_x[2]]), βbnds_x[1]])

            W_series[ttt] = W_series[ttt-1] + (1 - exp(-ξ*Δτ))*(W_bar-W_series[ttt-1]) + exp(-ξ*Δτ)*((endo_wealth*(W_series[ttt-1]*r_series[ttt-1] + Ereturn - drift_W))*Δτ + sqrt(Δτ)*η_r*rshocks[ttt-1] + sqrt(Δτ)*η_β*βshocks[ttt-1])
            W_series[ttt] = maximum([minimum([W_series[ttt], exp(Wbnds_x[2])]), exp(Wbnds_x[1])])

        end
        
        # store state series for histograms
        @views r_series_all[:,nnn] .= r_series[1:end- T_burn]
        @views β_series_all[:,nnn] .= β_series[1:end- T_burn]
        @views W_series_all[:,nnn] .= W_series[1:end- T_burn]
    end

    return SimulationQ(vec(r_series_all), vec(β_series_all), vec(W_series_all))

end


function normalize1(vec::Vector{R}, limit::I)  where {R <: Real, I <: Integer}
    vecn = vec[2:limit] .- vec[1]
    return vecn
end


function normalize2(mat::Matrix{R}) where {R <: Real}
    matn = mat[:,2] - mat[:,1]
    return matn
end


function calcIR(prms::params{R,I}, pitp, sitp, aEitp, sigmaitp, initvec, shocks, sim, wshock, sample_dummy, ir0)  where {R <: Real, I <: Integer}

    @unpack Wbnds, Wbnds_x, σ_r_unc, T_sim, N_sample, T_burn, αvec, θvec, W_bar, y_idxs, n_idxs, τ_y, n_τ, T_ir, W_center,
            mu, r_bar, β_bar, κ_r, κ_β, σ_r, σ_β, Δτ, θ0vec, θvec, τ, μ_W, τ_max = prms

    if sample_dummy == 1
        N_ir = N_sample
    else
        N_ir = 1
    end

    r_mat        = zeros(T_ir,N_ir)
    β_mat        = zeros(T_ir,N_ir)
    W_mat        = zeros(T_ir,N_ir)
    y_mat        = zeros(n_τ,T_ir,N_ir)
    aE_mat       = zeros(Int(τ_max), Int(τ_max), T_ir, N_ir)
    
    # sampled from simulation 
    rng = MersenneTwister(1234);
    start_idx_vec = rand(rng, 1:T_sim - T_burn - T_ir,N_ir)

    Threads.@threads for sss = 1:N_ir

        if N_ir == 1
            initstate = initvec
            rshocks = zeros(T_ir-1)
            βshocks = zeros(T_ir-1)
            Wshocks = zeros(T_ir-1) .+ wshock
        else 
            start_idx = start_idx_vec[sss]
            initstate = [sim.r[start_idx],sim.β[start_idx], maximum([minimum([sim.W[start_idx], exp(Wbnds_x[2]) - 0.05*W_center ]), exp(Wbnds_x[1]) + 0.05*W_center ]) ]
            rshocks = zeros(T_ir-1)
            βshocks = zeros(T_ir-1)
            Wshocks = zeros(T_ir-1)  .+ wshock           
        end

        rshocks[1] = rshocks[1] + shocks[1]
        βshocks[1] = βshocks[1] + shocks[2]
        Wshocks[1] = Wshocks[1] + shocks[3]

        r_mat[:,sss], β_mat[:,sss] , W_mat[:,sss], y_mat[:,:,sss], aE_mat[:,:,:,sss] = calcSeries(prms, initstate, rshocks, βshocks, Wshocks, pitp, sitp, aEitp, sigmaitp, false) 

    end

    r_series  = dropdims(mean(r_mat,dims=2),dims=2)
    β_series  = dropdims(mean(β_mat,dims=2),dims=2)
    W_series  = dropdims(mean(W_mat,dims=2),dims=2)
    y_series  = dropdims(mean(y_mat,dims=3),dims=3)
    aE_series = dropdims(mean(aE_mat,dims=4),dims=4)

    if ir0 != []
        r_series  = r_series[1]      .+ (r_series - ir0.r) 
        β_series  = β_series[1]      .+ (β_series - ir0.β)
        W_series  = W_series[1]      .+ (W_series - ir0.W) 
        y_series  = y_series[:,1]    .+ (y_series - ir0.y)
        aE_series = aE_series[:,:,1] .+ (aE_series - ir0.aE)
    end 
    yy_series = y_series[y_idxs,:] 

    tp_series = y_series - ((r_series' .- r_bar).*((1 .- exp.(-κ_r.*τ))./κ_r) .+ r_bar.*τ)./τ
    yr5f5_series = (10*y_series[y_idxs[10],:] - 5*y_series[y_idxs[5],:])./5
    tp5f5_series = yr5f5_series - ((r_series .- r_bar).*((exp.(-κ_r.*5) .- exp.(-κ_r.*10))./κ_r) .+ r_bar.*5)./5
    x_series = -(αvec.*y_series.*τ) .+ θ0vec .+ adj_fun.(W_series',Ref(prms)).*θvec.*β_series'
    fwd_series, ret_series, mret_series = calcRetFwd(prms, y_series, yy_series)
    CPss, CPshock   = CPDecompIR(prms, y_series, fwd_series, ret_series)
    CPssE, CPshockE = CPEDecompIR(prms, fwd_series, aE_series)

    return Impulseresponse(r_series, β_series, W_series, y_series, aE_series, fwd_series, ret_series, mret_series, [CPss, CPshock], [CPssE, CPshockE], tp_series, tp5f5_series, x_series)

end


function calcSeries(prms::params{R,I}, initstate::Vector{R}, rshocks::Vector{R}, βshocks::Vector{R}, Wshocks::Vector{R}, pitp, sitp, aEitp, sigmaitp, simdummy::Bool) where {R <: Real, I <: Integer}

    @unpack  ξ, W_bar, r_bar, β_bar, κ_r, κ_β,  θ0vec, θvec, αvec,  σ_r, σ_β, τ_max, rbnds, rbnds_x, βbnds, βbnds_x, Wbnds, Wbnds_x, Δτ, n_τ, τ, endo_wealth, n_dims, calc_aE_sim = prms

    n_t       = size(rshocks,1)+1      # length of sample  
    T         = Δτ*(n_t-1)  

    # preallocation
    r_series  = zeros(n_t)             # path of interest rates
    β_series  = zeros(n_t)             # path of demand
    W_series  = zeros(n_t)             # path of wealth
    y_series  = zeros(n_τ, n_t)
    aE_series = zeros(Int(τ_max), Int(τ_max), n_t)
    sigma_series  = zeros(Int(τ_max), n_t)
    sim_p    = zeros(n_τ,n_t)
    
    # starting point
    r_series[1]  = initstate[1]
    β_series[1]  = initstate[2]
    W_series[1]  = initstate[3]

    if n_dims == 3

        for τττ = 1:n_τ
            sim_p[τττ,1] = pitp[τττ](r_series[1],β_series[1],log(W_series[1]))
        end
        for τττ = 1:Int(τ_max)
            if simdummy == false || calc_aE_sim == true
                for tttttt = 1:Int(τ_max)
                    aE_series[τττ,tttttt,1] = aEitp[τττ,tttttt](r_series[1],β_series[1],log(W_series[1]))
                end
            end
            sigma_series[τττ,1] = sigmaitp[τττ](r_series[1],β_series[1],log(W_series[1]))
        end
    else
        for τττ = 1:n_τ
            sim_p[τττ,1] = pitp[τττ](r_series[1],β_series[1])
        end
        for τττ = 1:Int(τ_max)
            if simdummy == false || calc_aE_sim == true
                for tttttt = 1:Int(τ_max)
                    aE_series[τττ,tttttt,1] = aEitp[τττ,tttttt](r_series[1],β_series[1])
                end
            end
            sigma_series[τττ,1] = sigmaitp[τττ](r_series[1],β_series[1])
        end
    end

    for ttt = 2:n_t
        
        if n_dims == 3 
            drifts = sitp(r_series[ttt-1],β_series[ttt-1],log(W_series[ttt-1]))
        else 
            drifts = sitp(r_series[ttt-1],β_series[ttt-1])
        end

        Ereturn = drifts.myvec[4] # dot(coeffs_Er, polys)
        η_r     = drifts.myvec[5] # dot(coeffs_ηr, polys)
        η_β     = drifts.myvec[6] # dot(coeffs_ηβ, polys)

        r_series[ttt] = r_series[ttt-1] + κ_r*(r_bar - r_series[ttt-1])*Δτ + sqrt(Δτ)*σ_r*rshocks[ttt-1]
        r_series[ttt] = maximum([minimum([r_series[ttt], rbnds_x[2]]), rbnds_x[1]])
        β_series[ttt] = β_series[ttt-1] + κ_β*(β_bar - β_series[ttt-1])*Δτ + sqrt(Δτ)*σ_β*βshocks[ttt-1]
        β_series[ttt] = maximum([minimum([β_series[ttt], βbnds_x[2]]), βbnds_x[1]])

        W_series[ttt] = W_series[ttt-1] + (1 - exp(-ξ*Δτ))*(W_bar-W_series[ttt-1]) + exp(-ξ*Δτ)*((endo_wealth*(W_series[ttt-1]*r_series[ttt-1] + Ereturn))*Δτ + sqrt(Δτ)*η_r*rshocks[ttt-1] + sqrt(Δτ)*η_β*βshocks[ttt-1]) + Wshocks[ttt-1]
        W_series[ttt] = maximum([minimum([W_series[ttt], exp(Wbnds_x[2])]), exp(Wbnds_x[1])])

        if n_dims == 3
            for τττ = 1:n_τ
                sim_p[τττ,ttt] = pitp[τττ](r_series[ttt],β_series[ttt],log(W_series[ttt]))
            end
            for τττ = 1:Int(τ_max)
                if simdummy == false || calc_aE_sim == true
                    for tttttt = 1:Int(τ_max)
                        aE_series[τττ,tttttt,ttt] = aEitp[τττ,tttttt](r_series[ttt],β_series[ttt],log(W_series[ttt]))
                    end
                end
                sigma_series[τττ,ttt] = sigmaitp[τττ](r_series[ttt],β_series[ttt],log(W_series[ttt]))
            end
        else
            for τττ = 1:n_τ
                sim_p[τττ,ttt] = pitp[τττ](r_series[ttt],β_series[ttt])
            end
            if simdummy == false || calc_aE_sim == true
                for tttttt = 1:Int(τ_max)
                    for τττ = 1:Int(τ_max)
                        aE_series[τττ,tttttt,ttt] = aEitp[τττ,tttttt](r_series[ttt],β_series[ttt])
                    end
                end
            end
            for τττ = 1:Int(τ_max)
                sigma_series[τττ,ttt] = sigmaitp[τττ](r_series[ttt],β_series[ttt])
            end
        end

    end 

    y_series[2:end,:] = -sim_p[2:end,:]./τ[2:end]

    return r_series, β_series, W_series, y_series, aE_series, sigma_series

end


function calcIRR(prms::params{R,I}, pitp, sitp, aEitp, initvec, wshock, ir0, r_bar_shift)  where {R <: Real, I <: Integer}

    @unpack ξ, σ_r_unc, T_sim, T_burn, αvec, θvec, W_bar, y_idxs, n_idxs, τ_y, n_τ, T_irr, convergeCrit,
         mu, r_bar, β_bar, κ_r, κ_β, σ_r, σ_β, Δτ, θ0vec, θvec, τ, n_τ, μ_W, endo_wealth, Wbnds, Wbnds_x, rbnds, rbnds_x, βbnds, βbnds_x, τ_max, n_dims = prms

    r_series  = zeros(T_irr)                            # path of interest rates
    β_series  = zeros(T_irr)                            # path of demand shocks
    W_series  = zeros(T_irr)                            # path of wealth
    y_series = zeros(n_τ, T_irr)                        # path of yields
    aE_series = zeros(Int(τ_max), Int(τ_max), T_irr)    # path of expected returns

    r_series[1] = initvec[1]
    β_series[1] = initvec[2]
    W_series[1] = initvec[3]

    if n_dims == 3
        for τττ = 2:n_τ
            y_series[τττ,1] = -pitp[τττ](r_series[1],β_series[1],log(W_series[1]))/τ[τττ] + r_bar_shift[1]
        end
        for τττ = 1:Int(τ_max)
            for tttttt = 1:Int(τ_max)
                aE_series[τττ,tttttt,1] = aEitp[τττ,tttttt](r_series[1],β_series[1],log(W_series[1])) + r_bar_shift[1]
            end
        end
    else
        for τττ = 2:n_τ
            y_series[τττ,1] = -pitp[τττ](r_series[1],β_series[1])/τ[τττ] + r_bar_shift[1]
        end
        for τττ = 1:Int(τ_max)
            for tttttt = 1:Int(τ_max)
                aE_series[τττ,tttttt,1] = aEitp[τττ,tttttt](r_series[1],β_series[1]) + r_bar_shift[1]
            end
        end
    end

    for ttt = 2:T_irr

        r_series[ttt] = r_series[ttt-1] + κ_r*(r_bar - r_series[ttt-1])*Δτ - 0*(r_bar_shift[ttt]-r_bar_shift[ttt-1])
        r_series[ttt] = maximum([minimum([r_series[ttt], rbnds_x[2]]), rbnds_x[1]])
        β_series[ttt] = β_series[ttt-1] + κ_β*(β_bar - β_series[ttt-1])*Δτ
        β_series[ttt] = maximum([minimum([β_series[ttt], βbnds_x[2]]), βbnds_x[1]])
        
        wnxt_temp = W_series[ttt-1] + (1.0 - exp(-ξ*Δτ))*(W_bar - W_series[ttt-1])    # initial condition for next period wealth

        pnxt      = zeros(n_τ)                                                        # next period prices
        p_temp    = -y_series[:,ttt-1].*τ                                             # current period prices
        
        dif = 1.0
        while dif > convergeCrit
          
            # prices tomorrow
            for τττ = 2:n_τ
                pnxt[τττ] = pitp[τττ]([r_series[ttt],β_series[ttt],log(wnxt_temp)][1:n_dims]...) - r_bar_shift[ttt]*τ[τττ]
            end
            
            wnxt_new =  W_series[ttt-1] + (1.0 - exp(-ξ*Δτ))*(W_bar - W_series[ttt-1]) + exp(-ξ*Δτ)*((endo_wealth == 1)*(W_series[ttt-1]*r_series[ttt-1]*Δτ + 
                       Δτ*(pnxt[1:end-1] - p_temp[2:end] .- (r_series[ttt-1] .+ r_bar_shift[ttt-1])*Δτ)'*(αvec[2:end].*(p_temp[2:end] + τ[2:end].*r_bar_shift[ttt-1]) + θ0vec[2:end] + adj_fun(W_series[ttt-1],prms).*β_series[ttt-1].*θvec[2:end]))) + wshock
            
            if mu[3]>1
                if wnxt_new<exp(Wbnds_x[1])
                    wnxt_new = exp(Wbnds_x[1])
                end
                if wnxt_new>exp(Wbnds_x[2])
                    wnxt_new = exp(Wbnds_x[2])
                end
                dif = maximum(abs.(wnxt_new - wnxt_temp))
                wnxt_temp = deepcopy(wnxt_new)
            else
                dif = 0.0
            end    
        
        end

        W_series[ttt] = wnxt_temp
        y_series[2:end,ttt] = -pnxt[2:end]./τ[2:end]

        for τττ = 1:Int(τ_max)
            for tttttt = 1:Int(τ_max)
                aE_series[τττ,tttttt,ttt] = aEitp[τττ,tttttt]([r_series[ttt],β_series[ttt],log(W_series[ttt])][1:n_dims]...) + r_bar_shift[ttt]
            end
        end

    end

    if ir0 != []
        r_series  = r_series[1]   .+ (r_series - ir0.r[1:T_irr]) .+ r_bar_shift
        β_series  = β_series[1]   .+ (β_series - ir0.β[1:T_irr])
        W_series  = W_series[1]   .+ (W_series - ir0.W[1:T_irr])
        y_series  = y_series[:,1] .+ (y_series - ir0.y[:,1:T_irr])
        aE_series = aE_series[:,:,1] .+ (aE_series - ir0.aE[:,:,1:T_irr])
    end
    yy_series =y_series[y_idxs,:]
    tp_series = y_series - ((r_series' .- r_bar .- r_bar_shift').*((1 .- exp.(-κ_r.*τ))./κ_r) .+ (r_bar .+ r_bar_shift').*τ)./τ
    yr5f5_series = (10*y_series[y_idxs[10],:] - 5*y_series[y_idxs[5],:])./5
    tp5f5_series = yr5f5_series - ((r_series .- r_bar .- r_bar_shift).*((exp.(-κ_r.*5) .- exp.(-κ_r.*10))./κ_r) .+ (r_bar .+ r_bar_shift).*5)./5
    fwd_series, ret_series, mret_series = calcRetFwd(prms, y_series, yy_series)
    CPss, CPshock   = [], []
    CPssE, CPshockE = [], []
    x_series = []

    return Impulseresponse(r_series, β_series, W_series, y_series, aE_series, fwd_series, ret_series, mret_series, [CPss, CPshock], [CPssE, CPshockE], tp_series, tp5f5_series, x_series)

end


function calcIRM(prms::params{R,I}, pitp, sitp, aEitp, pitp_vec, sitp_vec, aEitp_vec, r_shiftvec, β_shiftvec, initvec, shocks, sim, wshock, sample_dummy, ir0)  where {R <: Real, I <: Integer}

    # unpack relevant parameters
    @unpack ξ, n_τ, mu, κ_r, γ, σ2_r, α, αvec, θvec, θ0vec, endo_wealth, T_sim, T_burn, N_sample, n_t,
             W_bar, W_center, r_bar, δ_al, δ_th, τ_max, θ0, θ, n_idxs, y_idxs, n_q, μ_W, convergeCrit, n_dims,
             σ_r, σ_r_unc, σ_β, σ_β_unc, σ2_β, κ_β, β_bar, rbnds, rbnds_x, βbnds, βbnds_x, Wbnds, Wbnds_x, τ, Δτ, T_ir, Δt, T_paste = prms
    
    ## now sample for transitions
    if sample_dummy == 1
        N_ir = N_sample
        # sampled from simulation 
        rng = MersenneTwister(1234);
        start_idx_vec = rand(rng, 1:T_sim - T_burn - T_ir,N_ir)
    else
        N_ir = 1
    end
    rshocks = zeros(T_ir-1,N_ir)
    βshocks = zeros(T_ir-1,N_ir)
    Wshocks = zeros(T_ir-1,N_ir) .+ wshock

    r_mat        = zeros(T_ir,N_ir)
    β_mat        = zeros(T_ir,N_ir)
    W_mat        = zeros(T_ir,N_ir)
    y_mat        = zeros(n_τ,T_ir,N_ir)
    aE_mat       = zeros(Int(τ_max), Int(τ_max), T_ir, N_ir)

    
    Threads.@threads for sss = 1:N_ir

        if N_ir == 1
            initstate = initvec
        else 
            start_idx = start_idx_vec[sss]
            initstate = [sim.r[start_idx],sim.β[start_idx], maximum([minimum([sim.W[start_idx], exp(Wbnds_x[2]) - 0.05*W_center ]), exp(Wbnds_x[1]) + 0.05*W_center ]) ]
        end
    
        # next state
        rnxt      = initstate[1] + κ_r*(r_bar-initstate[1])*Δτ 
        βnxt      = initstate[2] + κ_β*(β_bar-initstate[2])*Δτ 

        wnxt_temp = initstate[3] + (1.0-exp(-ξ*Δτ))*(W_bar - initstate[3]) 
        

        pnxt      = zeros(n_τ)
        p_temp    = zeros(n_τ)
        
        # prices today, use normal price interpolation
        for ttt = 2:n_τ
            p_temp[ttt] = pitp[ttt]([initstate[1],initstate[2],log(initstate[3])][1:n_dims]...)
        end    
                
        if sss == 1
        @printf "Find  transition.\n"
        end

        dif = 1.0
        while dif > 0.1*convergeCrit
          
            # prices tomorrow
            if n_dims > 2 
                for ttt = 2:n_τ
                    pnxt[ttt] = pitp_vec[1][ttt](rnxt,βnxt,log(wnxt_temp))
                end
            else
                for ttt = 2:n_τ
                    pnxt[ttt] = pitp_vec[1][ttt](rnxt,βnxt)
                end
            end 

            wnxt_new =  initstate[3] + (1.0-exp(-ξ*Δτ))*(W_bar - initstate[3])*Δτ + exp(-ξ*Δτ)*((endo_wealth == 1)*(initstate[3]*initstate[1]*Δτ + 
           Δτ*(pnxt[1:end-1] - p_temp[2:end] .- initstate[1]*Δτ)'*(αvec[2:end].*p_temp[2:end] + θ0vec[2:end] + adj_fun(initstate[3],prms).*initstate[2].*θvec[2:end])))
            
            if mu[3]>1
                if wnxt_new<exp(Wbnds_x[1])
                    wnxt_new = exp(Wbnds_x[1])
                end
                if wnxt_new>exp(Wbnds_x[2])
                    wnxt_new = exp(Wbnds_x[2])
                end
                dif = maximum(abs.(wnxt_new - wnxt_temp))
                wnxt_temp = deepcopy(wnxt_new)
            else
                dif = 0.0
            end    
        
        end

        # prices tomorrow
        if n_dims > 2 
            for ttt = 2:n_τ
                pnxt[ttt] = pitp_vec[1][ttt](rnxt,βnxt,log(wnxt_temp))
            end
        else
            for ttt = 2:n_τ
                pnxt[ttt] = pitp_vec[1][ttt](rnxt,βnxt)
            end
        end 

        # create impulse response
        r_mat[1,sss] = initstate[1]
        β_mat[1,sss] = initstate[2]
        W_mat[1,sss] = initstate[3]
        y_mat[2:end,1,sss] = -p_temp[2:end]./τ[2:end]

        r_mat[2,sss] = rnxt
        β_mat[2,sss] = βnxt
        W_mat[2,sss] = wnxt_temp
        y_mat[2:end,2,sss] = -pnxt[2:end]./τ[2:end]

        for tttttt = 1:Int(τ_max)
            for τττ = 1:Int(τ_max)
                aE_mat[τττ,tttttt,1,sss] = aEitp[τττ,tttttt]([r_mat[1,sss],β_mat[1,sss],log(W_mat[1,sss])][1:n_dims]...)
                aE_mat[τττ,tttttt,2,sss] = aEitp_vec[1][τττ,tttttt]([r_mat[2,sss],β_mat[2,sss],log(W_mat[2,sss])][1:n_dims]...)
            end
        end

        for ttt = 3:T_ir
            
            if (ttt-2 <= T_paste)
                drifts = sitp_vec[ttt-2]([r_mat[ttt-1,sss],β_mat[ttt-1,sss],log(W_mat[ttt-1,sss])][1:n_dims]...)
                r_current = r_mat[ttt-1,sss] + r_shiftvec[(ttt-3)*n_t + 1]
            else 
                drifts = sitp([r_mat[ttt-1,sss],β_mat[ttt-1,sss],log(W_mat[ttt-1,sss])][1:n_dims]...)
                r_current = r_mat[ttt-1,sss]
            end
            Ereturn = drifts.myvec[4] 
            η_r     = drifts.myvec[5] 
            η_β     = drifts.myvec[6] 

            r_mat[ttt,sss] = r_mat[ttt-1,sss] + κ_r*(r_bar - r_mat[ttt-1,sss])*Δτ + sqrt(Δτ)*σ_r*rshocks[ttt-1,sss]
            r_mat[ttt,sss] = maximum([minimum([r_mat[ttt,sss], rbnds_x[2]]), rbnds_x[1]])
            β_mat[ttt,sss] = β_mat[ttt-1,sss] + κ_β*(β_bar - β_mat[ttt-1,sss])*Δτ + sqrt(Δτ)*σ_β*βshocks[ttt-1,sss]
            β_mat[ttt,sss] = maximum([minimum([β_mat[ttt,sss], βbnds_x[2]]), βbnds_x[1]])

            W_mat[ttt,sss] = W_mat[ttt-1,sss] +  (1-exp(-ξ*Δτ))*(W_bar - W_mat[ttt-1,sss]) + exp(-ξ*Δτ)*((endo_wealth*(W_mat[ttt-1,sss]*r_current + Ereturn))*Δτ  + sqrt(Δτ)*η_r*rshocks[ttt-1,sss] + sqrt(Δτ)*η_β*βshocks[ttt-1,sss]) + Wshocks[ttt-1,sss]
            W_mat[ttt,sss] = maximum([minimum([W_mat[ttt,sss], exp(Wbnds_x[2])]), exp(Wbnds_x[1])])
            
            for τττ = 2:n_τ
                if (ttt-1 <= T_paste)
                    tmp = pitp_vec[ttt-1][τττ]([r_mat[ttt,sss],β_mat[ttt,sss],log(W_mat[ttt,sss])][1:n_dims]...)
                else
                    tmp = pitp[τττ]([r_mat[ttt,sss],β_mat[ttt,sss],log(W_mat[ttt,sss])][1:n_dims]...)
                end
                y_mat[τττ,ttt,sss] = -tmp/τ[τττ]
            end
            for tttttt = 1:Int(τ_max)
                for τττ = 1:Int(τ_max)
                    if (ttt-1 <= T_paste)
                        aE_mat[τττ,tttttt,ttt,sss] = aEitp_vec[ttt-1][τττ,tttttt]([r_mat[ttt,sss],β_mat[ttt,sss],log(W_mat[ttt,sss])][1:n_dims]...)
                    else
                        aE_mat[τττ,tttttt,ttt,sss] = aEitp[τττ,tttttt]([r_mat[ttt,sss],β_mat[ttt,sss],log(W_mat[ttt,sss])][1:n_dims]...)
                    end
                end
            end
            

        end

    end        

    r_series  = dropdims(mean(r_mat,dims=2),dims=2)
    β_series  = dropdims(mean(β_mat,dims=2),dims=2)
    W_series  = dropdims(mean(W_mat,dims=2),dims=2)
    y_series  = dropdims(mean(y_mat,dims=3),dims=3)
    aE_series = dropdims(mean(aE_mat,dims=4),dims=4)

    if ir0 != []
        r_series = r_series[1]   .+ (r_series - ir0.r) 
        β_series = β_series[1]   .+ (β_series - ir0.β)
        W_series = W_series[1]   .+ (W_series - ir0.W) 
        y_series = y_series[:,1] .+ (y_series - ir0.y)  
        aE_series = aE_series[:,:,1] .+ (aE_series - ir0.aE)
    end 

    yy_series = y_series[y_idxs,:] 
    fwd_series, ret_series, mret_series = calcRetFwd(prms, y_series, yy_series)
    CPss, CPshock = CPDecompIR(prms, yy_series, fwd_series, ret_series)
    CPssE, CPshockE = CPEDecompIR(prms, fwd_series, aE_series)
    
    # store shift in r
    r_series[2:T_paste+1] = r_series[2:T_paste+1] + r_shiftvec[1:n_t:T_paste*n_t]

    # store shift in habitat demand
    β_series[1] = sum((θ0vec).*Δτ,dims=1)[1]
    if size(β_shiftvec,2)<=size(β_series,1)-1
        β_series[2:size(β_shiftvec,2)+1] = sum((θ0vec .- β_shiftvec).*Δτ,dims=1)[1,:]
    else
        β_series[2:end] = sum((θ0vec - β_shiftvec[:,1:size(β_series,1)-1]).*Δτ,dims=1)[1,:]
    end

    tp_series = y_series - ((r_series' .- r_bar).*((1 .- exp.(-κ_r.*τ))./κ_r) .+ r_bar.*τ)./τ
    yr5f5_series = (10*y_series[y_idxs[10],:] - 5*y_series[y_idxs[5],:])./5
    tp5f5_series = yr5f5_series - ((r_series .- r_bar).*((exp.(-κ_r.*5) .- exp.(-κ_r.*10))./κ_r) .+ r_bar.*5)./5
    x_series = -(αvec.*y_series.*τ) .+ θ0vec .+ adj_fun.(W_series',Ref(prms)).*θvec.*β_series'

    return Impulseresponse(r_series, β_series, W_series, y_series, aE_series, fwd_series, ret_series, mret_series, [CPss, CPshock], [CPssE, CPshockE], tp_series, tp5f5_series, x_series)
    
end


function calcRetFwd(prms::params{R,I}, y_series::Matrix{R}, yy_series::Matrix{R})  where {R <: Real, I <: Integer}

    @unpack τ_y, Δτ_y, Δt_y, n_idxs, n_τ, Δτ, τ = prms

    n_t = size(y_series,2)

    ret_series = zeros(n_idxs-1,n_t)
    mret_series = zeros(n_τ-1,n_t)

    @views fwd_series = (τ_y[2:end].*yy_series[2:end,:] .- τ_y[1:end-1].*yy_series[1:end-1,:])./Δτ_y

    # annual return series
    @views ret_series[:,1:Δt_y]    .= NaN # can't calculate for first year of observations
    @views ret_series[:,Δt_y+1:end] = (-τ_y[1:end-1].*yy_series[1:end-1,Δt_y+1:end] .+ τ_y[2:end].*yy_series[2:end,1:end-Δt_y])./Δτ_y

    # period to period return series
    @views mret_series[:,1]        .= NaN
    @views mret_series[:,2:end]     = (-τ[1:end-1].*y_series[1:end-1,2:end] .+ τ[2:end].*y_series[2:end,1:end-1])./Δτ

    return fwd_series, ret_series, mret_series

end

function FBReg(prms::params{R,I}, yy_series::Matrix{R}, fwd_series::Matrix{R}, ret_series::Matrix{R} )  where {R <: Real, I <: Integer}

    @unpack n_idxs, Δt_y  = prms
    n_t = size(yy_series, 2)

    β_FB = zeros(n_idxs-1,1)
    for nnn = 1:n_idxs-1
        X = [ones(n_t-Δt_y,1) (fwd_series[nnn,1:end-Δt_y].-yy_series[1,1:end-Δt_y])]
        y = ret_series[nnn,Δt_y+1:end]-yy_series[1,1:end-Δt_y]
        β = (X'*X)\(X'*y)
        β_FB[nnn,1] = β[2,1]
    end

    return β_FB

end


function CSReg(prms::params{R,I}, yy_series::Matrix{R})  where {R <: Real, I <: Integer}

    @unpack Δt_y, n_idxs, y_idxs = prms

    n_t  = size(yy_series,2)-Δt_y

    β_CS = zeros(n_idxs-1,1)
    @simd for nnn = 2:n_idxs
        @views X = (1/(nnn-1)).*(yy_series[nnn, 1:n_t].-yy_series[1, 1:n_t])
        X = [ones(n_t,1) X]
        @views y = yy_series[nnn-1, Δt_y+1:end] .- yy_series[nnn, 1:n_t]
        β = (X'*X)\(X'*y)
        β_CS[nnn-1,1] = β[2,1]
    end

    return β_CS

end


function calcRegIRM(prms::params{R,I}, pitp, pitp_vec, pitp0_vec, initvec, sim, ppp)  where {R <: Real, I <: Integer}

    # unpack relevant parameters
    @unpack ξ, n_τ, mu, κ_r, γ, σ2_r, α, αvec, θvec, θ0vec, endo_wealth, T_sim, T_burn, N_sample, n_t,
             W_bar, W_center, r_bar, δ_al, δ_th, τ_max, θ0, θ, n_idxs, y_idxs, n_q, μ_W, convergeCrit, n_dims,
             σ_r, σ_r_unc, σ_β, σ_β_unc, σ2_β, κ_β, β_bar, rbnds, rbnds_x, βbnds, βbnds_x, Wbnds, Wbnds_x, τ, Δτ, T_ir, Δt, T_paste = prms
    
    ## now sample for transitions
    N_ir = N_sample
    rng = MersenneTwister(1234);
    start_idx_vec = rand(rng, 1:T_sim - T_burn - T_ir,N_ir)

    r_mat    = zeros(2,N_ir)
    β_mat    = zeros(2,N_ir)
    W_mat    = zeros(2,N_ir)
    y_mat    = zeros(n_τ,2,N_ir)
    f_mat    = zeros(n_idxs,2,N_ir)
    tp_mat   = zeros(2,N_ir)
    dur_mat  = zeros(2,N_ir)
    ig_mat   = zeros(2,N_ir)

    r0_mat   = zeros(2,N_ir)
    β0_mat   = zeros(2,N_ir)
    W0_mat   = zeros(2,N_ir)
    y0_mat   = zeros(n_τ,2,N_ir)
    f0_mat   = zeros(n_idxs,2,N_ir)
    tp0_mat  = zeros(2,N_ir)
    dur0_mat = zeros(2,N_ir)
    ig0_mat  = zeros(2,N_ir)

    # find transition
    for sss = 1:N_ir

        if N_ir == 1
            initstate = initvec
        else 
            start_idx = start_idx_vec[sss]
            initstate = [sim.r[start_idx],sim.β[start_idx], maximum([minimum([sim.W[start_idx], exp(Wbnds_x[2]) - 0.05*W_center ]), exp(Wbnds_x[1]) + 0.05*W_center ]) ]
        end
    
        # next state
        rnxt      = initstate[1] + κ_r*(r_bar-initstate[1])*Δτ 
        βnxt      = initstate[2] + κ_β*(β_bar-initstate[2])*Δτ 

        wnxt_temp  = initstate[3] + (1.0-exp(-ξ*Δτ))*(W_bar - initstate[3])
        wnxt0_temp = initstate[3] + (1.0-exp(-ξ*Δτ))*(W_bar - initstate[3])

        pnxt      = zeros(n_τ)
        pnxt0     = zeros(n_τ)
        p_temp    = zeros(n_τ)
        
        # prices today, use normal price interpolation
        if n_dims > 2 
            for ttt = 2:n_τ
               p_temp[ttt] = pitp[ttt](initstate[1],initstate[2],log(initstate[3]))
            end    
        else 
            for ttt = 2:n_τ
               p_temp[ttt] = pitp[ttt](initstate[1],initstate[2])
            end    
        end
                
        if sss == 1
        @printf "Find  transition.\n"
        end

        # find transition in presence of shock
        dif = 1.0
        while dif > convergeCrit
          
            # prices tomorrow
            if n_dims > 2 
                for ttt = 2:n_τ
                    pnxt[ttt] = pitp_vec[1][ttt](rnxt,βnxt,log(wnxt_temp))
                end
            else
                for ttt = 2:n_τ
                    pnxt[ttt] = pitp_vec[1][ttt](rnxt,βnxt)
                end
            end 

            wnxt_new =  initstate[3] + (1.0-exp(-ξ*Δτ))*(W_bar - initstate[3])*Δτ + exp(-ξ*Δτ)*((endo_wealth == 1)*(initstate[3]*initstate[1]*Δτ + 
                        Δτ*(pnxt[1:end-1] - p_temp[2:end] .- initstate[1]*Δτ)'*(αvec[2:end].*p_temp[2:end] + θ0vec[2:end] + adj_fun(initstate[3],prms)*initstate[2].*θvec[2:end])))
            
            if mu[3]>1
                if wnxt_new<exp(Wbnds_x[1])
                    wnxt_new = exp(Wbnds_x[1])
                end
                if wnxt_new>exp(Wbnds_x[2])
                    wnxt_new = exp(Wbnds_x[2])
                end
                dif = maximum(abs.(wnxt_new - wnxt_temp))
                wnxt_temp = deepcopy(wnxt_new)
            else
                dif = 0.0
            end    
        
        end

        # Store state variables and yields

        r_mat[1,sss] = initstate[1]
        β_mat[1,sss] = initstate[2]
        W_mat[1,sss] = initstate[3]
        y_mat[2:end,1,sss] = -p_temp[2:end]./τ[2:end]

        r_mat[2,sss] = rnxt
        β_mat[2,sss] = βnxt
        W_mat[2,sss] = wnxt_temp
        y_mat[2:end,2,sss] = -pnxt[2:end]./τ[2:end]

        # find transition in absence of shock
        dif = 1.0
        while dif > convergeCrit
          
            # prices tomorrow
            if n_dims > 2
                for ttt = 2:n_τ
                    pnxt0[ttt] = pitp0_vec[1][ttt](rnxt,βnxt,log(wnxt0_temp))
                end
            else 
                for ttt = 2:n_τ
                    pnxt0[ttt] = pitp0_vec[1][ttt](rnxt,βnxt)
                end
            end 

            wnxt_new =  initstate[3] + (1.0-exp(-ξ*Δτ))*(W_bar - initstate[3])*Δτ + exp(-ξ*Δτ)*((endo_wealth == 1)*(initstate[3]*initstate[1]*Δτ + 
                        Δτ*(pnxt0[1:end-1] - p_temp[2:end] .- initstate[1]*Δτ)'*(αvec[2:end].*p_temp[2:end] + θ0vec[2:end] + adj_fun(initstate[3], prms)*initstate[2].*θvec[2:end])))
            
            if mu[3]>1
                if wnxt_new<exp(Wbnds_x[1])
                    wnxt_new = exp(Wbnds_x[1])
                end
                if wnxt_new>exp(Wbnds_x[2])
                    wnxt_new = exp(Wbnds_x[2])
                end
                dif = maximum(abs.(wnxt_new - wnxt0_temp))
                wnxt0_temp = deepcopy(wnxt_new)
            else
                dif = 0.0
            end    
        
        end

        # Store state variables and yields

        r0_mat[1,sss] = initstate[1]
        β0_mat[1,sss] = initstate[2]
        W0_mat[1,sss] = initstate[3]
        y0_mat[2:end,1,sss] = -p_temp[2:end]./τ[2:end]

        r0_mat[2,sss] = rnxt
        β0_mat[2,sss] = βnxt
        W0_mat[2,sss] = wnxt0_temp
        y0_mat[2:end,2,sss] = -pnxt0[2:end]./τ[2:end]

        # Store other important variables

        for ttt = 1:2

            # Store forward rates
            f_mat[1,ttt,sss] = y_mat[y_idxs[1],ttt,sss]
            f0_mat[1,ttt,sss] = y0_mat[y_idxs[1],ttt,sss]
            for τττ = 2:Int(τ_max)
                f_mat[τττ,ttt,sss] = τττ*y_mat[y_idxs[τττ],ttt,sss] - (τττ-1)*y_mat[y_idxs[τττ-1],ttt,sss]
                f0_mat[τττ,ttt,sss] = τττ*y0_mat[y_idxs[τττ],ttt,sss] - (τττ-1)*y0_mat[y_idxs[τττ-1],ttt,sss]
            end

            # Store term premium

            yr5f5_temp = (10*y_mat[y_idxs[10],ttt,sss] - 5*y_mat[y_idxs[5],ttt,sss])./5
            tp_mat[ttt,sss] = yr5f5_temp - ((r_mat[ttt,sss] .- r_bar).*((exp.(-κ_r.*5) .- exp.(-κ_r.*10))./κ_r) .+ r_bar.*5)./5
            yr5f5_temp = (10*y0_mat[y_idxs[10],ttt,sss] - 5*y0_mat[y_idxs[5],ttt,sss])./5
            tp0_mat[ttt,sss] = yr5f5_temp - ((r0_mat[ttt,sss] .- r_bar).*((exp.(-κ_r.*5) .- exp.(-κ_r.*10))./κ_r) .+ r_bar.*5)./5

            # Store IG and duration

            X_temp = -(αvec.*y_mat[:,ttt,sss].*τ) + θ0vec + adj_fun(W_mat[ttt,sss],prms)*θvec.*β_mat[ttt,sss]
            A_temp = zeros(size(X_temp,1))
            A_temp[X_temp.>=0] = X_temp[X_temp.>=0]

            A_temp  = (sum(A_temp.*Δτ) + (W_mat[ttt,sss] - sum(X_temp.*Δτ)).*(W_mat[ttt,sss] - sum(X_temp.*Δτ) .> 0))
            X1_temp = X_temp[y_idxs[1]+1:end]

            ig_mat[ttt,sss]  = ((W_mat[ttt,sss] - sum(X1_temp.*Δτ))/A_temp)
            dur_mat[ttt,sss] = (sum(X_temp.*τ.*Δτ)./W_mat[ttt,sss])

            X_temp = -(αvec.*y0_mat[:,ttt,sss].*τ) + θ0vec +adj_fun(W0_mat[ttt,sss],prms)*θvec.*β0_mat[ttt,sss]
            A_temp = zeros(size(X_temp,1))
            A_temp[X_temp.>=0] = X_temp[X_temp.>=0]

            A_temp  = (sum(A_temp.*Δτ) + (W0_mat[ttt,sss] - sum(X_temp.*Δτ)).*(W0_mat[ttt,sss] - sum(X_temp.*Δτ) .> 0))
            X1_temp = X_temp[y_idxs[1]+1:end]

            ig0_mat[ttt,sss]  = ((W0_mat[ttt,sss] - sum(X1_temp.*Δτ))/A_temp)
            dur0_mat[ttt,sss] = (sum(X_temp.*τ.*Δτ)./W0_mat[ttt,sss])

        end

    end

    # Compute deltas
    dy1 = (y_mat[y_idxs[1],2,:] - y0_mat[y_idxs[1],2,:]).*100
    df  = (f_mat[:,2,:] - f0_mat[:,2,:]).*100

    # Compute other regressors
    xtp  = tp0_mat[2,:].*100
    xdur = dur0_mat[2,:]
    xig  = ig0_mat[2,:]
    yreg = df[20,:]
    yreg2 = yreg[xdur.>0]
    
    # Regression against 5y5 term premium
    
    Xreg = [ones(size(xtp,1)) xtp]
    βreg = (Xreg'*Xreg)\(Xreg'*yreg)
    βtp5f5 = βreg[2]/mean(dy1)

    # Regression against log duration
    Xreg = [ones(size(xdur[xdur.>0],1)) log.(xdur[xdur.>0])]
    βreg = (Xreg'*Xreg)\(Xreg'*yreg2)
    βdur = βreg[2]/mean(dy1)

    # Regression against income gap
    Xreg = [ones(size(xig,1)) -xig]
    βreg = (Xreg'*Xreg)\(Xreg'*yreg)
    βig = βreg[2]/mean(dy1)

    # Regression against 5y5 term premium (restricted)
    Xreg = [ones(size(xtp[xdur.>0],1)) xtp[xdur.>0]]
    βreg = (Xreg'*Xreg)\(Xreg'*yreg2)
    βtp5f5_rest = βreg[2]/mean(dy1)

    # Regression against income gap
    Xreg = [ones(size(xig[xdur.>0],1)) -xig[xdur.>0]]
    βreg = (Xreg'*Xreg)\(Xreg'*yreg2)
    βig_restr = βreg[2]/mean(dy1)
    
    return StateDepReg(βtp5f5, βdur, βig, βtp5f5_rest, βig_restr)
    
end


function CPDecompIR(prms::params{R,I}, y_series::Matrix{R}, fwd_series::Matrix{R}, ret_series::Matrix{R})  where {R <: Real, I <: Integer}

    @unpack Δt_y, y_idxs = prms

    n_τ      = size(fwd_series,1)

    fwd_ss = fwd_series[:,1]
    Ey1_ss = y_series[1,1]
    lhs_ss = fwd_ss .- Ey1_ss
    rhs_ss = zeros(n_τ, n_τ)
    for ccc = 1:n_τ
        for rrr = 1:n_τ
            if ccc < rrr
                rhs_ss[rrr,ccc] = (ret_series[rrr-ccc+1,1] - Ey1_ss) - (ret_series[rrr-ccc,1] - Ey1_ss)
            elseif ccc == rrr
                rhs_ss[rrr,ccc] = ret_series[1,1] - Ey1_ss
            end
        end
    end

    fwd = fwd_series[:,2]
    Ey1 = y_series[1,2+Δt_y:Δt_y:2+n_τ*Δt_y]
    lhs = fwd - Ey1
    rhs = zeros(n_τ, n_τ)
    for ccc = 1:n_τ
        for rrr = 1:n_τ
            if ccc < rrr
                rhs[rrr,ccc] = (ret_series[rrr-ccc+1,2+Δt_y*(ccc-1)] - y_series[1,2+Δt_y*(ccc-1)]) - (ret_series[rrr-ccc,2+Δt_y*(ccc-1)] - y_series[1,2+Δt_y*(ccc-1)])
            elseif ccc == rrr
                rhs[rrr,ccc] = ret_series[1,2+Δt_y*(rrr-1)] - y_series[1,2+Δt_y*(rrr-1)]
            end
        end
    end

    CPMatrix_ss = [lhs_ss rhs_ss]
    CPMatrix    = [lhs rhs]

    return CPMatrix_ss, CPMatrix

end


function CPDecomp(prms::params{R,I}, y_series::Matrix{R}, fwd_series::Matrix{R}, ret_series::Matrix{R})  where {R <: Real, I <: Integer}

    @unpack Δt_y, y_idxs, τ_max = prms

    n_τ       = size(fwd_series,1)
    t_len     = size(ret_series,2)

    y_series  = y_series[y_idxs,:]
    CP_series = zeros(Int(τ_max-1),Int(τ_max), t_len-n_τ*Δt_y)

    for ttt = 1:t_len-n_τ*Δt_y

        fwd_ss = fwd_series[:,ttt]                      # current yearly forward rates
        Ey1_ss = y_series[1,ttt+Δt_y:Δt_y:ttt+n_τ*Δt_y] # 1y yields over next n_tau years
        lhs_ss = fwd_ss .- Ey1_ss                       # spread forwards - 1y yields 
        rhs_ss = zeros(n_τ, n_τ)                        
        for ccc = 1:n_τ
            for rrr = 1:n_τ
                if ccc < rrr
                    rhs_ss[rrr,ccc] = (ret_series[rrr-ccc+1,ttt+Δt_y*(ccc-1)] - y_series[1,ttt+Δt_y*(ccc-1)]) - (ret_series[rrr-ccc,ttt+Δt_y*(ccc-1)] - y_series[1,ttt+Δt_y*(ccc-1)])
                elseif ccc == rrr
                    rhs_ss[rrr,ccc] = ret_series[1,ttt+Δt_y*(rrr-1)] - y_series[1,ttt+Δt_y*(rrr-1)]
                end
            end
        end

        CP_series[:,:,ttt] = [lhs_ss rhs_ss]
    end
    return CP_series

end


function CPEDecompIR(prms::params{R,I}, fwd_series::Matrix{R}, ret_series::Array{R})  where {R <: Real, I <: Integer}

    @unpack Δt_y, y_idxs = prms

    n_τ      = size(fwd_series,1)

    fwd_ss = fwd_series[:,1]
    Ey1_ss = ret_series[1,2:end,1]
    lhs_ss = fwd_ss .- Ey1_ss
    rhs_ss = zeros(n_τ, n_τ)
    for ccc = 1:n_τ
        for rrr = 1:n_τ
            if ccc < rrr
                rhs_ss[rrr,ccc] = ret_series[rrr-ccc+2,ccc,1] - ret_series[rrr-ccc+1,ccc,1]
            elseif ccc == rrr
                rhs_ss[rrr,ccc] = ret_series[2,rrr,1] - ret_series[1,rrr,1]
            end
        end
    end

    fwd = fwd_series[:,2]
    Ey1 = ret_series[1,2:end,2]
    lhs = fwd - Ey1
    rhs = zeros(n_τ, n_τ)
    for ccc = 1:n_τ
        for rrr = 1:n_τ
            if ccc < rrr
                rhs[rrr,ccc] = ret_series[rrr-ccc+2,ccc,2] - ret_series[rrr-ccc+1,ccc,2]
            elseif ccc == rrr
                rhs[rrr,ccc] = ret_series[2,rrr,2] - ret_series[1,rrr,2]
            end
        end
    end

    CPMatrix_ss = [lhs_ss rhs_ss]
    CPMatrix    = [lhs rhs]

    return CPMatrix_ss, CPMatrix

end


function CPEDecompIRR(prms::params{R,I}, fwd_series::Matrix{R}, ret_series::Array{R})  where {R <: Real, I <: Integer}

    @unpack Δt_y, y_idxs = prms

    n_t      = size(fwd_series,2)
    carry_mat = zeros(5, n_t)

    for ttt = 1:n_t
        for ccc = 1:5
            carry_mat[ccc,ttt] = (ret_series[10+1-ccc,ccc,ttt] - ret_series[5+1-ccc,ccc,ttt])/5
        end
    end

    return carry_mat

end


function CPEDecomp(prms::params{R,I}, fwd_series::Matrix{R}, ret_series::Array{R})  where {R <: Real, I <: Integer}

    @unpack Δt_y, y_idxs, τ_max= prms

    n_τ      = size(fwd_series,1)
    t_len    = size(fwd_series,2)
    CP_series = zeros(Int(τ_max-1),Int(τ_max), t_len-1)

    for ttt = 1:t_len-1

    fwd_ss = fwd_series[:,ttt]
    Ey1_ss = ret_series[1,2:end,ttt]
    lhs_ss = fwd_ss .- Ey1_ss
    rhs_ss = zeros(n_τ, n_τ)
    for ccc = 1:n_τ
        for rrr = 1:n_τ
            if ccc < rrr
                rhs_ss[rrr,ccc] = ret_series[rrr-ccc+2,ccc,ttt] - ret_series[rrr-ccc+1,ccc,ttt]
            elseif ccc == rrr
                rhs_ss[rrr,ccc] = ret_series[2,rrr,ttt] - ret_series[1,rrr,ttt]
            end
        end
    end

    CP_series[:,:,ttt] = [lhs_ss rhs_ss]

    end
    
    return CP_series

end
