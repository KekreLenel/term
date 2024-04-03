function adj_fun(W, prms::params{R,I}) where {R<:Real, I<:Integer}
    @unpack adj_coeff, W_center = prms 
    return 1.0- 1.0/(1+adj_coeff*W/W_center)
end

function calcP(prms::params{R,I}) where {R<:Real, I<:Integer}

    # unpack parameters
    @unpack τ, n_τ, n_mc, n_grid, n_q, μ_r, μ_β, μ_W, rbnds, βbnds, Wbnds, αvec, n_dims, convergeCrit,
            θ0vec, θvec, γ, σ_r, σ_β, Δτ, Δt, endo_wealth, n_t, interp_prices, interp_drifts = prms # dr, dβ, dW 
    
    # number of threads available
    n_thrds   = Threads.nthreads()

    # quadrature calculation
    nodes, weights  = gausshermite(n_q)
    weights         = weights/sqrt(pi)
    nodes           = nodes*sqrt(2)
    rquad           = kron(ones(n_q),nodes)
    βquad           = kron(nodes, ones(n_q))
    nquad           = n_q^2
    quadshocks      = [rquad βquad]'
    weights         = kron(weights,weights)

    # preallocation
    P            = ones(n_τ, n_grid)                         # price array
    P_dt         = ones(n_τ, n_grid)                         # price array in +dt
    P_r          = zeros(n_τ, n_grid)                        # price derivative r
    P_β          = zeros(n_τ, n_grid)                        # price derivative β
    P_W          = zeros(n_τ, n_grid)                        # price derivative W
    tmp_array    = [zeros(n_τ, n_grid) for _ in 1:n_thrds]   # tmp storage of simulated prices
    tmp_array_dt = [zeros(n_τ, n_grid) for _ in 1:n_thrds]   # tmp storage of simulated prices in +dt
    tmpE_array   = [zeros(n_τ, n_grid) for _ in 1:nquad]     # tmp storage of simulated price expectations
    EP           = zeros(n_τ, n_grid)                        # expected price array
    EdP_mc       = zeros(n_τ, n_grid)                        # expected log returns array

    demand_P_r  = zeros(n_grid)                             # sum (demand*P_r*dtau)
    demand_P_β  = zeros(n_grid)                             # sum (demand*P_β*dtau)
    demand_P_W  = zeros(n_grid)                             # sum (demand*P_W*dtau)
   
    drift_r_vec = zeros(n_grid)                             # risk-neutral r drift
    drift_β_vec = zeros(n_grid)                             # risk-neutral β drift    
    drift_W_vec = zeros(n_grid)                             # risk-neutral W drift    
    Ereturn_vec = zeros(n_grid)                             # expected portfolio returns    
    η_r_vec     = zeros(n_grid)                             # η_r loading of W on r
    η_β_vec     = zeros(n_grid)                             # η_β loading of W on β 

    drift_r_vec_new = zeros(n_grid)                         # updated of the above
    drift_β_vec_new = zeros(n_grid)                         # ''
    drift_W_vec_new = zeros(n_grid)                         # ''
    Ereturn_vec_new = zeros(n_grid)                         # ''
    η_r_vec_new     = zeros(n_grid)                         # ''
    η_β_vec_new     = zeros(n_grid)                         # ''

    # shocks for MC sims: antithetic sampling
    rng    = MersenneTwister(1234)
    shocks = [randn(rng, 2, (n_τ-1)*n_t) for _ in 1:div(n_mc,2)]
    shocks = vcat(shocks, -shocks)

    # grid over states
    r0     = range(rbnds..., μ_r)
    β0     = range(βbnds..., μ_β)
    logW0  = range(Wbnds..., μ_W)
    W0     = exp.(logW0)

    # collect all states as vector of static vectors
    transstates     = [SVector(r, β, W) for W in logW0 for β in β0 for r in r0 ]
    states          = [SVector(r, β, W) for W in W0 for β in β0 for r in r0 ]
    states_array    = reduce(hcat, states)'; # also store as array 

    # prepare vector of price interpolations
    itp       = interpolate(reshape(log.(P[1,:]), (μ_r, μ_β, μ_W)[1:n_dims]), interp_drifts) # interpolate
    etpf      = extrapolate(itp, Line())                                                      # choose whether to extrapolate
    sitp      = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)                      # scale interpolation
    
    itp       = interpolate(reshape(log.(P[1,:]), (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices) # interpolate
    etpf      = extrapolate(itp, Line())                                           # choose whether to extrapolate
    pitp      = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)                          # scale interpolation

    pitp      = Vector{typeof(pitp)}(undef,n_τ)                                    # prep vector of interpolations across maturities
    pitp_dt   = deepcopy(pitp)

    drift_vecs = [MyVec(SVector{6}(drift_r_vec[nnn], drift_β_vec[nnn], drift_W_vec[nnn], Ereturn_vec[nnn], η_r_vec[nnn], η_β_vec[nnn])) for nnn = 1:n_grid]

    t0 = time_ns()
    iter = 0
    dif = 1.0
    while dif > convergeCrit || iter < 10

        iter += 1
        t1 = time_ns()
        
        # reset storage arrays 
        tmp_array  .*= 0.0
        tmp_array_dt .*= 0.0
        tmpE_array .*= 0.0
       
        # drift and loading interpolation
        itp = interpolate(reshape(drift_vecs, (μ_r, μ_β, μ_W)[1:n_dims]), interp_drifts) 
        etpf = extrapolate(itp, Line())
        sitp = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)
        
        # monte carlo price simulation 
        Threads.@threads for mmm = 1:n_mc
            getSample!(tmp_array[Threads.threadid()], tmp_array_dt[Threads.threadid()], states, sitp, shocks[mmm], prms)
        end
    
        # average over MC sims
        P    .= sum(tmp_array)/n_mc
        P_dt .= sum(tmp_array_dt)/n_mc
        
        if any(isnan.(P))
            error("P is nan")
        end
        if any(P .<= sqrt(eps()))
            error("P is nan")
        end
        
        # interpolate over prices, get derivatives
        for ttt = 1:n_τ

            itp = interpolate(reshape(log.(P[ttt,:]), (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices)
            etpf      = extrapolate(itp, Line())
            pitp[ttt] = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)

            itp = interpolate(reshape(log.(P_dt[ttt,:]), (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices)
            etpf      = extrapolate(itp, Line())
            pitp_dt[ttt] = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)

            derivs = [Interpolations.gradient(pitp[ttt], transstates[nnn][1:n_dims]...) for nnn = 1:n_grid]
            derivs_temp = [derivs[nnn][1] for nnn = 1:n_grid]
            P_r[ttt,:] = P[ttt,:].*derivs_temp
            derivs_temp = [derivs[nnn][2] for nnn = 1:n_grid]
            P_β[ttt,:] = P[ttt,:].*derivs_temp
            if n_dims == 3
                derivs_temp = [derivs[nnn][3] for nnn = 1:n_grid]
                P_W[ttt,:] = P[ttt,:].*derivs_temp./states_array[:,3]
            else 
                P_W[ttt,:] .= 0.0
            end

        end
        
        # get price expectations for expected returns 
        for mmm = 1:nquad
            getSampleE!(tmpE_array[mmm], states, pitp_dt, 0.0, sitp, quadshocks[:,mmm], weights[mmm], prms)
        end
       
        EP .= sum(tmpE_array)

        if any(isnan.(EP))
            error("EP is nan")
        end

        # expected returns
        EdP_mc[2:end,:] = log.(EP[2:end,:]) -  log.(P[2:end,:]) 
       
        # demand 
        demand      = (αvec.*log.(P) .+ θ0vec .+ adj_fun.(states_array[:,3]',Ref(prms)).*states_array[:,2]'.*θvec)
        if any(isnan.(demand))
            error("demand is nan")
        end       

        # calculate relevant integrals
        demand_P_r  = sum(demand.*P_r./P, dims=1)'*Δτ
        demand_P_β  = sum(demand.*P_β./P, dims=1)'*Δτ
        demand_P_W  = sum(demand.*P_W./P, dims=1)'*Δτ
        if any(demand_P_W .> 0.9)
            println("Warning: demand_P_W > 0.9") # temporarily fine, but should not trigger at final iteration.
            demand_P_W[demand_P_W.>0.9] .= 0.9
        end 

        @views Ereturn_vec_new  = Δτ .*sum(demand.*(EdP_mc./Δt .- states_array[:,1]'), dims=1)'
        η_r_vec_new      = demand_P_r*σ_r./(1.0 .- demand_P_W)
        η_β_vec_new      = demand_P_β*σ_β./(1.0 .- demand_P_W)    

        @views drift_r_vec_new = γ./states_array[:,3].*σ_r.*(σ_r*demand_P_r + η_r_vec.*demand_P_W)
        @views drift_β_vec_new = γ./states_array[:,3].*σ_β.*(σ_β*demand_P_β + η_β_vec.*demand_P_W)
        @views drift_W_vec_new = γ./states_array[:,3].*(η_r_vec.*(σ_r*demand_P_r + η_r_vec.*demand_P_W) + η_β_vec.*(σ_β*demand_P_β + η_β_vec.*demand_P_W))
       
        if any(drift_r_vec_new .> 20.0)
            println("Warning: drift_r_vec_new > 20.0.") # temporarily fine, but should not trigger at final iteration.
            drift_r_vec_new[drift_r_vec_new.>20.0] .= 20.0
        end

        if any(isnan.(drift_r_vec_new))
            error("drift_r_vec_new is nan")
        end    
        if any(isnan.(drift_β_vec_new))
            error("drift_β_vec_new is nan")
        end    
        if any(isnan.(drift_W_vec_new))
            error("drift_W_vec_new is nan")
        end  
        if any(isnan.(Ereturn_vec_new))
            error("Ereturn_vec_new is nan")
        end
        if any(isnan.(η_r_vec_new))
            error("η_r_vec_new is nan")
        end
        if any(isnan.(η_β_vec_new))
            error("η_β_vec_new is nan")
        end     

        dif_vec = [ maximum(abs.(drift_r_vec - drift_r_vec_new)) , maximum(abs.(drift_β_vec - drift_β_vec_new)) , maximum(abs.(drift_W_vec - drift_W_vec_new)) , 
                   endo_wealth*maximum(abs.(Ereturn_vec - Ereturn_vec_new)) , endo_wealth*maximum(abs.(η_r_vec - η_r_vec_new)) , endo_wealth*maximum(abs.(η_β_vec - η_β_vec_new))]
        dif = sum(dif_vec)

        t2 = time_ns()
        
        if mod(iter,10) == 0
        @printf "Diff %0.3e %0.2e %0.2e %0.2e %0.2e %0.2e %0.2e %i %0.3fs\n" dif dif_vec... iter float((t2-t1)/1.0e9)
        end
        flush(stdout)

        # update drifts 
        drift_r_vec = drift_r_vec + 0.1*(drift_r_vec_new-drift_r_vec)
        drift_β_vec = drift_β_vec + 0.1*(drift_β_vec_new-drift_β_vec)
        drift_W_vec = drift_W_vec + 0.05*(drift_W_vec_new-drift_W_vec)
        Ereturn_vec = Ereturn_vec + endo_wealth*0.2*(Ereturn_vec_new-Ereturn_vec)
        η_r_vec     = η_r_vec     + endo_wealth*0.05*(η_r_vec_new-η_r_vec)
        η_β_vec     = η_β_vec     + endo_wealth*0.05*(η_β_vec_new-η_β_vec)

        drift_vecs = [MyVec(SVector{6}(drift_r_vec[nnn], drift_β_vec[nnn], drift_W_vec[nnn], Ereturn_vec[nnn], η_r_vec[nnn], η_β_vec[nnn])) for nnn = 1:n_grid]

        if iter == 750
            println("Warning: No convergence!") 
            dif = 0.0
        end

        
    end 

    t2 = time_ns()
    @printf "\nConvergence! Final diff: %0.3e %i\n" dif iter
    @printf "Time to convergence: %0.3fs\n\n" float((t2-t0)/1.0e9)
    flush(stdout)

    return pitp, sitp, P, states, drift_vecs, r0, β0, logW0, nquad, weights, quadshocks, states, shocks, drift_r_vec, drift_β_vec, drift_W_vec, Ereturn_vec, η_r_vec, η_β_vec

end

function calcPx(prms::params{R,I}, pitp_noshift, sitp, drift_vecs, nquad, weights, quadshocks, rshift, QEpurchase, W_avrg) where {R<:Real, I<:Integer}

    # unpack parameters
    @unpack τ, n_τ, n_mc, n_grid, μ_r, μ_β, μ_W, rbnds, βbnds, Wbnds, αvec, κ_m, n_dims, τ_max, 
            interp_prices, interp_drifts, θ0vec, θvec, γ, σ_r, σ_β, Δτ, Δt, n_t, T_paste, gdp, y_idxs = prms 

    r_shiftvec    = zeros(T_paste*n_t)
    if rshift != [] 
        r_shiftvec[1] = rshift*σ_r*sqrt(Δτ)
        for rrr = 2:T_paste*n_t
            r_shiftvec[rrr] = r_shiftvec[rrr-1] + κ_m*(0. .- r_shiftvec[rrr-1])*Δt
        end   
    end
    
    β_shiftvec = zeros(n_τ, T_paste)
    if QEpurchase != [] 
            
        β_shiftvec[2:y_idxs[30],1] = QEpurchase[:,1].*(W_avrg/(0.1*gdp*Δτ))
        for bbb = 2:T_paste
            β_shiftvec[2:end-1,bbb] = deepcopy(β_shiftvec[3:end,bbb-1])
            if bbb <= size(QEpurchase,2)
                β_shiftvec[2:y_idxs[30],bbb] = β_shiftvec[2:y_idxs[30],bbb] + QEpurchase[:,bbb].*(1.0/(0.29*gdp*Δτ))# (W_avrg/(0.19*gdp*Δτ))
            end
        end

    end

    # shocks for MC sims
    rng    = MersenneTwister(1234)
    shocks = [randn(rng, 2, (n_τ-1)*n_t) for _ in 1:div(n_mc,2)]
    shocks = vcat(shocks, -shocks)

    # number of threads available
    n_thrds   = Threads.nthreads()

    # grid over states
    r0     = range(rbnds..., μ_r)
    β0     = range(βbnds..., μ_β)
    logW0  = range(Wbnds..., μ_W)
    W0     = exp.(logW0)

    # collect all states as vector of static vectors
    transstates     = [SVector(r, β, W) for W in logW0 for β in β0 for r in r0 ]
    states          = [SVector(r, β, W) for W in W0 for β in β0 for r in r0 ]
    states_array    = reduce(hcat, states)'; # also store as array 

    # preallocation
    P              = ones(n_τ, n_grid)                         # price array
    P_dt           = ones(n_τ, n_grid)                         # price array
    P_r            = zeros(n_τ, n_grid)                        # price derivative r
    P_β            = zeros(n_τ, n_grid)                        # price derivative β
    P_W            = zeros(n_τ, n_grid)                        # price derivative W
    aE_ret         = zeros(Int(τ_max), n_grid, Int(τ_max))

    tmp_array      = [zeros(n_τ, n_grid) for _ in 1:n_thrds]                    # tmp storage of simulated prices
    tmp_array_dt   = [zeros(n_τ, n_grid) for _ in 1:n_thrds]                    # tmp storage of simulated prices
    tmpE_array     = [zeros(n_τ, n_grid) for _ in 1:nquad]                      # tmp storage of simulated price expectations
    tmpaE_array    = [zeros(Int(τ_max), n_grid, Int(τ_max)) for _ in 1:n_thrds] # tmp storage of simulated price expectations with annual frequency

    EP          = zeros(n_τ, n_grid)                        # expected price array
    EdP_mc      = zeros(n_τ, n_grid)                        # expected log returns array

    demand_P_r  = zeros(n_grid)                             # sum (demand*P_r*dtau)
    demand_P_β  = zeros(n_grid)                             # sum (demand*P_β*dtau)
    demand_P_W  = zeros(n_grid)                             # sum (demand*P_W*dtau)
   
    drift_r_vec = zeros(n_grid)                             # risk-neutral r drift
    drift_β_vec = zeros(n_grid)                             # risk-neutral β drift    
    drift_W_vec = zeros(n_grid)                             # risk-neutral W drift    
    Ereturn_vec = zeros(n_grid)                             # expected portfolio returns    
    η_r_vec     = zeros(n_grid)                             # η_r loading of W on r
    η_β_vec     = zeros(n_grid)                             # η_β loading of W on β 

    
    itp       = interpolate(reshape(log.(P[1,:]), (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices) # interpolate
    etpf      = extrapolate(itp, Line())                                                     # choose whether to extrapolate
    pitp      = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)                     # scale interpolation
    pitp      = Vector{typeof(pitp)}(undef,n_τ)                                              # prep vector of interpolations across maturities
    pitp_dt   = deepcopy(pitp)

    itp       = interpolate(reshape(aE_ret[1,:,1], (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices) # interpolate
    etpf      = extrapolate(itp, Line())                                                      # choose whether to extrapolate
    aEitp     = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)                      # scale interpolation
    aEitp     = Array{typeof(aEitp)}(undef, Int(τ_max), Int(τ_max))                        # prep vector of interpolations across maturities
        
    drift_mat = [drift_vecs for _ = 1:T_paste]
    sitp_vec  = [sitp for _ = 1:T_paste]
    pitp_vec  = [pitp for _ = 1:T_paste]
    aEitp_vec = [aEitp for _ = 1:T_paste]

    # calculate prices and drift recursively
    t0 = time_ns()
    for rrr = T_paste:-1:1

        t1 = time_ns()
        
        tmp_array  .*= 0.0
        tmp_array_dt .*= 0.0
        tmpE_array .*= 0.0
        tmpaE_array .*= 0.0

        Threads.@threads for mmm = 1:n_mc
            getSample_shifted!(tmp_array[Threads.threadid()], tmp_array_dt[Threads.threadid()], states, sitp, sitp_vec, rrr,  shocks[mmm], r_shiftvec, prms)
        end

        P    .= sum(tmp_array)/n_mc
        P_dt .= sum(tmp_array_dt)/n_mc

        if any(isnan.(P))
            error("P is nan")
        end
        if any(isnan.(EP))
            error("EP is nan")
        end
        if any(P .<= sqrt(eps()))
            error("P is nan")
        end

        for ttt = 1:n_τ
            itp = interpolate(reshape(log.(P[ttt,:]), (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices)
            etpf      = extrapolate(itp, Line())
            pitp[ttt] = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)

            itp = interpolate(reshape(log.(P_dt[ttt,:]), (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices)
            etpf      = extrapolate(itp, Line())
            pitp_dt[ttt] = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)

            derivs = [Interpolations.gradient(pitp[ttt], transstates[nnn][1:n_dims]...) for nnn = 1:n_grid]
            derivs_temp = [derivs[nnn][1] for nnn = 1:n_grid]
            P_r[ttt,:] = P[ttt,:].*derivs_temp
            derivs_temp = [derivs[nnn][2] for nnn = 1:n_grid]
            P_β[ttt,:] = P[ttt,:].*derivs_temp
            if n_dims == 3
                derivs_temp = [derivs[nnn][3] for nnn = 1:n_grid]
                P_W[ttt,:] = P[ttt,:].*derivs_temp./states_array[:,3]
            else 
                P_W[ttt,:] = zeros(n_grid)
            end
        end

        for mmm = 1:nquad
            getSampleE!(tmpE_array[mmm], states, pitp_dt, r_shiftvec[(rrr-1)*n_t + 1], sitp_vec[rrr], quadshocks[:,mmm], weights[mmm], prms)
        end
        
        EP .= sum(tmpE_array)

        # calculate expected returns
        EdP_mc[2:end,:] = log.(EP[2:end,:]) -  log.(P[2:end,:]) 
            
        demand      = (αvec.*log.(P) .+ θ0vec .+ adj_fun.(states_array[:,3]',Ref(prms)).*states_array[:,2]'.*θvec .- β_shiftvec[:,rrr])
        if any(isnan.(demand))
            error("demand is nan")
        end       

        # calculate relevant integrals
        demand_P_r  = sum(demand.*P_r./P*Δτ, dims=1)'
        demand_P_W  = sum(demand.*P_W./P*Δτ, dims=1)'
        demand_P_β  = sum(demand.*P_β./P*Δτ, dims=1)'
        demand_P_W[demand_P_W.>0.9] .= 0.9

        @views Ereturn_vec  = Δτ .*sum(demand.*(EdP_mc./Δt .- (states_array[:,1]' .+ r_shiftvec[(rrr-1)*n_t+1] ) ), dims=1)'
        η_r_vec      = demand_P_r*σ_r./(1.0 .- demand_P_W)
        η_β_vec      = demand_P_β*σ_β./(1.0 .- demand_P_W)    

        @views drift_r_vec = γ./states_array[:,3].*σ_r.*(σ_r*demand_P_r + η_r_vec.*demand_P_W)
        @views drift_β_vec = γ./states_array[:,3].*σ_β.*(σ_β*demand_P_β + η_β_vec.*demand_P_W)
        @views drift_W_vec = γ./states_array[:,3].*(η_r_vec.*(σ_r*demand_P_r + η_r_vec.*demand_P_W) + η_β_vec.*(σ_β*demand_P_β + η_β_vec.*demand_P_W))
        drift_W_vec[drift_W_vec.>20.0] .= 20.0

        if any(isnan.(drift_r_vec))
            error("drift_r_vec_new is nan")
        end    
        if any(isnan.(drift_β_vec))
            error("drift_β_vec_new is nan")
        end    
        if any(isnan.(drift_W_vec))
            error("drift_W_vec_new is nan")
        end  
        if any(isnan.(Ereturn_vec))
            error("Ereturn_vec_new is nan")
        end
        if any(isnan.(η_r_vec))
            error("η_r_vec_new is nan")
        end
        if any(isnan.(η_β_vec))
            error("η_β_vec_new is nan")
        end

        if rrr > 1 
            drift_mat[rrr-1] = [MyVec(SVector{6}(drift_r_vec[nnn], drift_β_vec[nnn], drift_W_vec[nnn], Ereturn_vec[nnn], η_r_vec[nnn], η_β_vec[nnn])) for nnn = 1:n_grid]

            itp = interpolate(reshape(drift_mat[rrr-1], (μ_r, μ_β, μ_W)[1:n_dims]), interp_drifts ) 
            etpf = extrapolate(itp, Line())
            sitp_vec[rrr-1] = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)
        end

        # convert price array into interpolating array
        for ttt = 1:n_τ
            itp       = interpolate(reshape(log.(P[ttt,:]), (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices)
            etpf      = extrapolate(itp, Line())
            pitp[ttt] = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)
        end
        pitp_vec[rrr] = deepcopy(pitp)

        
        # compute annual expected prices up to 30 years ahead
        Threads.@threads for mmm = 1:n_mc
            getSampleECarry_shifted!(tmpaE_array[Threads.threadid()], states, pitp_noshift, pitp_vec, sitp, sitp_vec, shocks[mmm], rrr, r_shiftvec, prms)
        end
        aEp = sum(tmpaE_array)/n_mc

        # annual expected returns, today and up to 30 years ahead
        aE_ret = zeros(Int(τ_max), n_grid, Int(τ_max))
        @views aE_ret[1,:,1] .= -log.(P[y_idxs[1],:])
        @views aE_ret[1,:,2:end] .= -aEp[1,:,1:end-1]
        @views aE_ret[2:end,:,1] .= aEp[1:end-1,:,1] .- log.(P[y_idxs[2:end],:])
        @views aE_ret[2:end,:,2:end] .= aEp[1:end-1,:,2:end] .- aEp[2:end,:,1:end-1]
        
        # convert annual expected returns array into interpolating array
        @simd for tttttt = 1:Int(τ_max)
            @simd for ttt = 1:Int(τ_max)

                if n_dims > 2 
                    @views itp  = interpolate(reshape(aE_ret[ttt,:,tttttt], (μ_r, μ_β, μ_W)), interp_prices)
                    etpf = extrapolate(itp, Line())
                    aEitp[ttt,tttttt] = Interpolations.scale(etpf, (r0, β0, logW0))
                else
                    @views itp  = interpolate(reshape(aE_ret[ttt,:,tttttt], (μ_r, μ_β)), interp_prices)
                    etpf = extrapolate(itp, Line())
                    aEitp[ttt,tttttt] = Interpolations.scale(etpf, (r0, β0))
                end
            
            end
    
        end
        
        aEitp_vec[rrr] = deepcopy(aEitp)
        

        t2 = time_ns()
        if mod((T_paste-rrr+1),10) == 0
        @printf "Recursion %i/%i %0.3fs\n" (T_paste-rrr+1) T_paste float((t2-t1)/1.0e9)
        end
        flush(stdout)

    end
    t2 = time_ns()
    @printf "\nDone. Total time %0.3fs\n" float((t2-t0)/1.0e9)
    flush(stdout)
    return pitp_vec, sitp_vec, aEitp_vec, r_shiftvec, β_shiftvec


end

function getSample!(tmp_array::AbstractArray{R}, tmp_array_dt::AbstractArray{R}, states::AbstractVector{SVector{3, R}}, sitp, shocks, prms::params{R,I}) where {R<:Real, I<:Integer}
    
    @unpack n_grid = prms
    @simd for ggg = 1:n_grid 
        getP!(@view(tmp_array[:,ggg]), @view(tmp_array_dt[:,ggg]), states[ggg], sitp, shocks, prms)
    end

end

function getSample_shifted!(tmp_array::AbstractArray{R}, tmp_array_dt::AbstractArray{R}, states::AbstractVector{SVector{3, R}}, sitp, sitp_vec, rrr, shocks, rshift ,prms::params{R,I}) where {R<:Real, I<:Integer}
    
    @unpack n_grid = prms
    for ggg = 1:n_grid 
        getP_shifted!(@view(tmp_array[:,ggg]), @view(tmp_array_dt[:,ggg]), states[ggg], sitp, sitp_vec, rrr, shocks, rshift, prms)
    end

end


function getSampleE!(tmp_array::AbstractArray{R}, states::AbstractVector{SVector{3, R}}, pitp, rshift, sitp, shocks, weight, prms::params{R,I}) where {R<:Real, I<:Integer}

    @unpack n_grid = prms
    for ggg = 1:n_grid 
        getEP!(@view(tmp_array[:,ggg]), states[ggg], pitp, rshift, sitp, shocks, weight, prms)
    end

end

function getSampleECarry!(tmp_array::AbstractArray{R}, states::AbstractVector{SVector{3, R}}, pitp, sitp, shocks, prms::params{R,I}) where {R<:Real, I<:Integer}

    @unpack n_grid = prms
    for ggg = 1:n_grid 
        getEPCarry!(@view(tmp_array[:,ggg,:]), states[ggg], pitp, sitp, shocks, prms)
    end

end

function getSampleE2Carry!(tmp_array::AbstractArray{R}, tmp2_array::AbstractArray{R}, states::AbstractVector{SVector{3, R}}, pitp, sitp, shocks, prms::params{R,I}) where {R<:Real, I<:Integer}

    @unpack n_grid = prms
    for ggg = 1:n_grid 
        getEP2Carry!(@view(tmp_array[:,ggg]), @view(tmp2_array[:,ggg]), states[ggg], pitp, sitp, shocks, prms)
    end

end

function getSampleECarry_shifted!(tmp_array::AbstractArray{R}, states::AbstractVector{SVector{3, R}}, pitp, pitp_vec, sitp, sitp_vec, shocks, rrr, rshift, prms::params{R,I}) where {R<:Real, I<:Integer}

    @unpack n_grid = prms
    for ggg = 1:n_grid 
        getEPCarry_shifted!(@view(tmp_array[:,ggg,:]), states[ggg], pitp, pitp_vec, sitp, sitp_vec, shocks, rrr, rshift, prms)
    end

end


function getP_shifted!(P_vec::AbstractArray{R}, P_dt_vec::AbstractArray{R}, states::AbstractVector{R}, sitp, sitp_vec, rrr, shocks::AbstractArray{R}, rshift::AbstractVector{R}, prms::params{R,I} )  where {R<:Real, I<:Integer}

    N = size(P_vec,1)

    @unpack σ_r, κ_r, ξ, r_bar, σ_β, κ_β, β_bar, Δt, n_t, W_bar, endo_wealth, Wbnds_x, βbnds_x, rbnds_x, T_paste, n_dims  = prms

    r = states[1]
    β = states[2]
    W = states[3]

    # starting period 
    P  = 1.0
    P_vec[1]  += P
    P_dt_vec[1]  += P

    # store first-period shift 
    first_shift = rshift[(rrr-1)*n_t + 1]

    for i in 2:(N-1)*n_t

        store = (mod(i-1,n_t) == 0)
        
        if i-1 <= (T_paste-(rrr-1))*n_t    
            current_r = (r+rshift[(rrr-1)*n_t + i-1])
        else
            current_r = r
        end

        P  *= exp(-current_r*Δt)

        # P_vec stores at monthly maturities, 
        if store 
            P_vec[Int((i-1)/n_t)+1] += P
            P_dt_vec[Int((i-1)/n_t)+1] += P*exp((first_shift + r)*Δt)
        end 

        if i-1 <= (T_paste-rrr)*n_t
            if n_dims == 3
                drifts = sitp_vec[Int(floor((i-2)/n_t + rrr))](r,β,log(W))
            else
                drifts = sitp_vec[Int(floor((i-2)/n_t + rrr))](r,β)
            end
        else
            if n_dims == 3
                drifts = sitp(r,β,log(W))
            else
                drifts = sitp(r,β)
            end
        end

        drift_r = drifts.myvec[1]
        drift_β = drifts.myvec[2]
        drift_W = drifts.myvec[3]
        Ereturn = drifts.myvec[4]
        η_r     = drifts.myvec[5]
        η_β     = drifts.myvec[6]
        
        W += (1-exp(-ξ*Δt))*(W_bar-W) + exp(-ξ*Δt)*((endo_wealth*W*current_r + Ereturn - drift_W)*Δt  + sqrt(Δt)*(η_r*shocks[1,i] + η_β*shocks[2,i]))
        W = W*(exp(Wbnds_x[1]) <= W <= exp(Wbnds_x[2])) + exp(Wbnds_x[1])*(exp(Wbnds_x[1]) > W) + exp(Wbnds_x[2])*(exp(Wbnds_x[2]) < W)
        
        r += (κ_r*(r_bar-r)-drift_r)*Δt + sqrt(Δt)*σ_r*shocks[1,i]
        r = r*(rbnds_x[1] <= r <= rbnds_x[2]) + rbnds_x[1]*(rbnds_x[1] > r) + rbnds_x[2]*(rbnds_x[2] < r)

        β += (κ_β*(β_bar-β)-drift_β)*Δt + sqrt(Δt)*σ_β*shocks[2,i]
        β = β*(βbnds_x[1] <= β <= βbnds_x[2]) + βbnds_x[1]*(βbnds_x[1] > β) + βbnds_x[2]*(βbnds_x[2] < β)

    end

    if (N-1)*n_t <= (T_paste-(rrr-1))*n_t    
        current_r = (r+rshift[(rrr-1)*n_t + (N-1)*n_t])
    else
        current_r = r
    end

    P *= exp(-current_r*Δt)
    P_vec[N] += P
    P_dt_vec[N] += P*exp((first_shift + r)*Δt)

end

function getP!( P_vec::AbstractArray{R}, P_dt_vec::AbstractArray{R}, states::AbstractVector{R}, sitp,  shocks::AbstractArray{R}, prms::params{R,I} )  where {R<:Real, I<:Integer}

    @unpack n_dims, ξ, σ_r, κ_r, r_bar, σ_β, κ_β, β_bar, Δt, n_t, W_bar, endo_wealth, Wbnds_x, βbnds_x, rbnds_x = prms

    N = size(P_vec,1)

    r = states[1]
    β = states[2]
    W = states[3]

    # starting period 
    P  = 1.0
    P_vec[1]  += P
    P_dt_vec[1] += P

    for i in 2:(N-1)*n_t # n_t = 10, N = 4 , 30x

        store = (mod(i-1,n_t) == 0)

        # P_dt_vec stores monthly maturities less dt
        if store 
            @inbounds P_dt_vec[Int((i-1)/n_t)+1] += P
        end
        
        P *= exp(-r*Δt)

        # P_vec stores at monthly maturities, 
        if store 
            @inbounds P_vec[Int((i-1)/n_t)+1] += P
        end 
        
        if n_dims > 2
            drifts = sitp(r,β,log(W))
        else
            drifts = sitp(r,β)
        end

        drift_r = drifts.myvec[1] 
        drift_β = drifts.myvec[2] 
        drift_W = drifts.myvec[3] 
        Ereturn = drifts.myvec[4] 
        η_r     = drifts.myvec[5] 
        η_β     = drifts.myvec[6] 
        
        W += (1-exp(-ξ*Δt))*(W_bar-W) + exp(-ξ*Δt)*((endo_wealth*W*r + Ereturn - drift_W)*Δt  + sqrt(Δt)*(η_r*shocks[1,i] + η_β*shocks[2,i]))
        W = W*(exp(Wbnds_x[1]) <= W <= exp(Wbnds_x[2])) + exp(Wbnds_x[1])*(exp(Wbnds_x[1]) > W) + exp(Wbnds_x[2])*(exp(Wbnds_x[2]) < W)
        
        r += (κ_r*(r_bar-r)-drift_r)*Δt + sqrt(Δt)*σ_r*shocks[1,i]
        r = r*(rbnds_x[1] <= r <= rbnds_x[2]) + rbnds_x[1]*(rbnds_x[1] > r) + rbnds_x[2]*(rbnds_x[2] < r)

        β += (κ_β*(β_bar-β)-drift_β)*Δt + sqrt(Δt)*σ_β*shocks[2,i]
        β = β*(βbnds_x[1] <= β <= βbnds_x[2]) + βbnds_x[1]*(βbnds_x[1] > β) + βbnds_x[2]*(βbnds_x[2] < β)


    end

    P_dt_vec[N] += P
    P *= exp(-r*Δt)
    P_vec[N] += P

end

# calculated expected prices dt periods from today
function getEP!(P_vec::AbstractArray{R}, states::AbstractVector{R}, pitp, rshift, sitp, shocks::AbstractArray{R}, weight, prms::params{R,I} ) where {R<:Real, I<:Integer}

    @unpack n_dims, ξ, σ_r, κ_r, r_bar, σ_β, κ_β, β_bar, Δt, W_bar, endo_wealth, Wbnds_x, βbnds_x, rbnds_x, n_τ, n_t  = prms

    r = states[1]
    β = states[2]
    W = states[3]

    # physical measure transition one period ahead (no drift adjustment) !!!
    if n_dims > 2
       drifts = sitp(r,β,log(W))
    else
       drifts = sitp(r,β)
    end

    Ereturn = drifts.myvec[4]
    η_r     = drifts.myvec[5] 
    η_β     = drifts.myvec[6] 
    
    W += (1-exp(-ξ*Δt))*(W_bar-W) + exp(-ξ*Δt)*((endo_wealth*W*(r+rshift) + Ereturn)*Δt  + sqrt(Δt)*(η_r*shocks[1,1] + η_β*shocks[2,1]))
    W = W*(exp(Wbnds_x[1]) <= W <= exp(Wbnds_x[2])) + exp(Wbnds_x[1])*(exp(Wbnds_x[1]) > W) + exp(Wbnds_x[2])*(exp(Wbnds_x[2]) < W)

    r += (κ_r*(r_bar-r))*Δt + sqrt(Δt)*σ_r*shocks[1,1]
    β += (κ_β*(β_bar-β))*Δt + sqrt(Δt)*σ_β*shocks[2,1]
    r = r*(rbnds_x[1] <= r <= rbnds_x[2]) + rbnds_x[1]*(rbnds_x[1] > r) + rbnds_x[2]*(rbnds_x[2] < r)
    β = β*(βbnds_x[1] <= β <= βbnds_x[2]) + βbnds_x[1]*(βbnds_x[1] > β) + βbnds_x[2]*(βbnds_x[2] < β)

    @simd for τττ = 1:n_τ
        if n_dims > 2
            @inbounds P_vec[τττ] += weight*exp.(pitp[τττ](r,β,log(W)))
        else
            @inbounds P_vec[τττ] += weight*exp.(pitp[τττ](r,β))
        end
    end


end

function getEPCarry!(p_vec::AbstractArray{R}, states::AbstractVector{R}, pitp, sitp, shocks::AbstractArray{R}, prms::params{R,I} )  where {R<:Real, I<:Integer}

    @unpack n_dims, ξ, σ_r, κ_r, r_bar, σ_β, κ_β, β_bar, Δt, n_t, W_bar, endo_wealth, Wbnds_x, βbnds_x, rbnds_x, n_idxs, y_idxs, n_τ, calccarry, f_τ = prms

    r = states[1]
    β = states[2]
    W = states[3]

    for i = 1:n_t*n_τ
        
        if i <= n_t*f_τ || calccarry == true

            if n_dims > 2
                drifts = sitp(r,β,log(W))
            else
                drifts = sitp(r,β)
            end

            Ereturn = drifts.myvec[4]
            η_r     = drifts.myvec[5]
            η_β     = drifts.myvec[6] 
            
            W += (1-exp(-ξ*Δt))*(W_bar-W) + exp(-ξ*Δt)*((endo_wealth*W*r + Ereturn)*Δt  + sqrt(Δt)*(η_r*shocks[1,i] + η_β*shocks[2,i]))
            W = W*(exp(Wbnds_x[1]) <= W <= exp(Wbnds_x[2])) + exp(Wbnds_x[1])*(exp(Wbnds_x[1]) > W) + exp(Wbnds_x[2])*(exp(Wbnds_x[2]) < W)

            r += (κ_r*(r_bar-r))*Δt + sqrt(Δt)*σ_r*shocks[1,i]
            β += (κ_β*(β_bar-β))*Δt + sqrt(Δt)*σ_β*shocks[2,i]
            r = r*(rbnds_x[1] <= r <= rbnds_x[2]) + rbnds_x[1]*(rbnds_x[1] > r) + rbnds_x[2]*(rbnds_x[2] < r)
            β = β*(βbnds_x[1] <= β <= βbnds_x[2]) + βbnds_x[1]*(βbnds_x[1] > β) + βbnds_x[2]*(βbnds_x[2] < β)

            if mod(i, n_t) == 0.0
                if n_dims > 2
                    for τττ = 1:n_τ
                        p_vec[τττ, Int(i/n_t)] += pitp[τττ](r,β,log(W))
                    end
                else
                    for τττ = 1:n_τ
                        p_vec[τττ, Int(i/n_t)] += pitp[τττ](r,β)
                    end
                end
            end

        end

    end

end

function getEP2Carry!(y_vec::AbstractArray{R}, y2_vec::AbstractArray{R}, states::AbstractVector{R}, pitp, sitp, shocks::AbstractArray{R}, prms::params{R,I} )  where {R<:Real, I<:Integer}

    @unpack ξ, n_dims, σ_r, κ_r, r_bar, σ_β, κ_β, β_bar, Δt, n_t, τ, τ_max, W_bar, endo_wealth, Wbnds_x, βbnds_x, rbnds_x, n_idxs, y_idxs, n_τ, calccarry, f_τ = prms

    r = states[1]
    β = states[2]
    W = states[3]

    for i = 1:n_t
        
        if n_dims > 2
            drifts = sitp(r,β,log(W))
        else
            drifts = sitp(r,β)
        end

        Ereturn = drifts.myvec[4]
        η_r     = drifts.myvec[5]
        η_β     = drifts.myvec[6] 
        
        W += (1-exp(-ξ*Δt))*(W_bar-W) + exp(-ξ*Δt)*((endo_wealth*W*r + Ereturn)*Δt  + sqrt(Δt)*(η_r*shocks[1,i] + η_β*shocks[2,i]))
        W = W*(exp(Wbnds_x[1]) <= W <= exp(Wbnds_x[2])) + exp(Wbnds_x[1])*(exp(Wbnds_x[1]) > W) + exp(Wbnds_x[2])*(exp(Wbnds_x[2]) < W)

        r += (κ_r*(r_bar-r))*Δt + sqrt(Δt)*σ_r*shocks[1,i]
        β += (κ_β*(β_bar-β))*Δt + sqrt(Δt)*σ_β*shocks[2,i]
        r = r*(rbnds_x[1] <= r <= rbnds_x[2]) + rbnds_x[1]*(rbnds_x[1] > r) + rbnds_x[2]*(rbnds_x[2] < r)
        β = β*(βbnds_x[1] <= β <= βbnds_x[2]) + βbnds_x[1]*(βbnds_x[1] > β) + βbnds_x[2]*(βbnds_x[2] < β)

    end

    if n_dims > 2
        for τττ = 1:Int(τ_max)
            y2_vec[τττ] += (-pitp[y_idxs[τττ]](r,β,log(W))/τ[τττ])^2
            y_vec[τττ] += -pitp[y_idxs[τττ]](r,β,log(W))/τ[τττ]
        end
    else
        for τττ = 1:Int(τ_max)
            y2_vec[τττ] += (-pitp[y_idxs[τττ]](r,β)/τ[τττ])^2
            y_vec[τττ] += -pitp[y_idxs[τττ]](r,β)/τ[τττ]
        end
    end

end

function getEPCarry_shifted!(p_vec::AbstractArray{R}, states::AbstractVector{R}, pitp, pitp_vec, sitp, sitp_vec, shocks::AbstractArray{R}, rrr, rshift::AbstractVector{R}, prms::params{R,I} )  where {R<:Real, I<:Integer}

    @unpack ξ, n_dims, σ_r, κ_r, r_bar, σ_β, κ_β, β_bar, Δt, n_t, W_bar, endo_wealth, Wbnds_x, βbnds_x, rbnds_x, n_idxs, y_idxs, n_τ, τ_max, T_paste, f_τ, calccarry_shocks = prms

    r = states[1]
    β = states[2]
    W = states[3]

    for i = 1:n_t*(n_τ-1)

        if rrr == 1 || i <= n_t*f_τ || calccarry_shocks == true 
        
            if i <= (T_paste-rrr+1)*n_t
                if n_dims == 3
                    drifts = sitp_vec[Int(floor((i-1)/n_t + rrr))](r,β,log(W))
                else
                    drifts = sitp_vec[Int(floor((i-1)/n_t + rrr))](r,β)
                end
                shift = rshift[(rrr-1)*n_t + i]
            else
                if n_dims == 3
                    drifts = sitp(r,β,log(W))
                else
                    drifts = sitp(r,β)
                end
                shift = 0.0
            end

            Ereturn = drifts.myvec[4]
            η_r     = drifts.myvec[5]
            η_β     = drifts.myvec[6] 

            W += (1-exp(-ξ*Δt))*(W_bar-W) + exp(-ξ*Δt)*((endo_wealth*W*(r + shift) + Ereturn)*Δt  + sqrt(Δt)*(η_r*shocks[1,i] + η_β*shocks[2,i]))
            W = W*(exp(Wbnds_x[1]) <= W <= exp(Wbnds_x[2])) + exp(Wbnds_x[1])*(exp(Wbnds_x[1]) > W) + exp(Wbnds_x[2])*(exp(Wbnds_x[2]) < W)

            r += (κ_r*(r_bar-r))*Δt + sqrt(Δt)*σ_r*shocks[1,i]
            β += (κ_β*(β_bar-β))*Δt + sqrt(Δt)*σ_β*shocks[2,i]
            r = r*(rbnds_x[1] <= r <= rbnds_x[2]) + rbnds_x[1]*(rbnds_x[1] > r) + rbnds_x[2]*(rbnds_x[2] < r)
            β = β*(βbnds_x[1] <= β <= βbnds_x[2]) + βbnds_x[1]*(βbnds_x[1] > β) + βbnds_x[2]*(βbnds_x[2] < β)

            if mod(i, n_t*f_τ) == 0.0
                if i <= (T_paste-rrr)*n_t
                    if n_dims > 2
                        for τττ = 1:Int(τ_max)
                            p_vec[τττ, Int(i/(n_t*f_τ))] += pitp_vec[Int(i/n_t + rrr)][y_idxs[τττ]](r,β,log(W))
                        end
                    else
                        for τττ = 1:Int(τ_max)
                            p_vec[τττ, Int(i/(n_t*f_τ))] += pitp_vec[Int(i/n_t + rrr)][y_idxs[τττ]](r,β)
                        end
                    end
                else
                    if n_dims > 2
                        for τττ = 1:Int(τ_max)
                            p_vec[τττ, Int(i/(n_t*f_τ))] += pitp[y_idxs[τττ]](r,β,log(W))
                        end
                    else
                        for τττ = 1:Int(τ_max)
                            p_vec[τττ, Int(i/(n_t*f_τ))] += pitp[y_idxs[τττ]](r,β)
                        end
                    end
                end
            end

        end

    end

end

function calcCarry(prms::params{R,I}, P, pitp, sitp, states, r0, β0, logW0)  where {R<:Real, I<:Integer}

    # unpack relevant parameters
    @unpack n_τ, mu, κ_r, γ, σ2_r, α, αvec, θvec, θ0vec, μ_r, n_mc, f_τ,
            W_bar, δ_al, δ_th, τ_max, θ0, θ, n_idxs, y_idxs, n_q, μ_W, μ_β, 
            σ_r, σ_β, σ2_β, κ_β, β_bar, rbnds, βbnds, Wbnds, τ, Δτ, n_grid, n_t, n_dims, interp_prices = prms

    # number of threads available
    n_thrds   = Threads.nthreads()
    
    # shocks for MC sims
    rng    = MersenneTwister(1234)
    shocks = [randn(rng, 2, n_t*n_τ) for _ in 1:div(n_mc,2)]
    shocks = vcat(shocks, -shocks)

    # preallocation
    tmpE_array = [zeros(n_τ, n_grid, n_τ) for _ in 1:n_thrds]  # tmp storage of simulated price expectations
    tmpy_array = [zeros(Int(τ_max), n_grid) for _ in 1:n_thrds]
    tmpy2_array = [zeros(Int(τ_max), n_grid) for _ in 1:n_thrds]

    # store expected prices at monthly frequency
    Threads.@threads for mmm = 1:n_mc
        getSampleECarry!(tmpE_array[Threads.threadid()], states, pitp, sitp, shocks[mmm], prms)
    end
    Ep = sum(tmpE_array)/n_mc

    # store expected squared prices next period
    Threads.@threads for mmm = 1:n_mc
        getSampleE2Carry!(tmpy_array[Threads.threadid()], tmpy2_array[Threads.threadid()], states, pitp, sitp, shocks[mmm], prms)
    end
    sigmay = (sum(tmpy2_array)/n_mc - (sum(tmpy_array)/n_mc).^2).^0.5

    # monthly expected returns, today and up to 30*12 months ahead
    E_ret = zeros(n_τ, n_grid, n_τ)
    E_ret[2:end,:,1] += Ep[1:end-1,:,1] - log.(P[2:end,:])
    E_ret[2:end,:,2:end] += Ep[1:end-1,:,2:end] - Ep[2:end,:,1:end-1]

    # annual expected returns, today and up to 30 years ahead
    aE_ret = zeros(Int(τ_max), n_grid, Int(τ_max))
    aE_ret[:,:,1] += Ep[y_idxs.-f_τ,:,y_idxs[1]-1] - log.(P[y_idxs,:])
    aE_ret[:,:,2:end] += Ep[y_idxs.-f_τ,:,y_idxs[2:end].-1] - Ep[y_idxs,:,y_idxs[1:end-1].-1]

    # carry returns and forward rates implied by carry returns
    carry_ret = E_ret[2:end,:, :] - E_ret[1:end-1, :, :]

    fwds = zeros(n_τ,n_grid)
    Threads.@threads for nnn = 1:n_grid
        tmp = zeros(n_τ-2,n_τ-2)
        for ttt = 1:(n_τ-2)
            for ttt2 = 1:(n_τ-2)
                if ttt<=ttt2
                    tmp[ttt2,ttt] = carry_ret[ttt2-ttt+2,nnn,ttt]
                end
            end
        end
        fwds[3:end,nnn] = sum(tmp, dims=2) + E_ret[2,nnn,2:n_τ-1]
        fwds[2,nnn] = E_ret[2,nnn,1]
    end

    # prepare matrix of interpolation coefficients for annual expected returns
    itp   = interpolate(reshape(aE_ret[1,:,1], (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices) # interpolate
    etpf  = extrapolate(itp, Line())                                                      # choose whether to extrapolate
    aEitp = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)                      # scale interpolation
    aEitp = Array{typeof(aEitp), 2}(undef, Int(τ_max), Int(τ_max))                        # prep vector of interpolations across maturities

    for ttt = 1:Int(τ_max)

        for tttttt = 1:Int(τ_max)

            itp  = interpolate(reshape(aE_ret[ttt,:,tttttt], (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices)
            etpf = extrapolate(itp, Line())
            aEitp[ttt,tttttt] = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims])
        
        end

    end

    # prepare matrix of interpolation coefficients for conditional volatility
    itp   = interpolate(reshape(sigmay[1,:], (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices) # interpolate
    etpf  = extrapolate(itp, Line())                                                      # choose whether to extrapolate
    sigmaitp = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims]...)                      # scale interpolation
    sigmaitp = Array{typeof(sigmaitp)}(undef, Int(τ_max))                        # prep vector of interpolations across maturities

    for ttt = 1:Int(τ_max)

        itp  = interpolate(reshape(sigmay[ttt,:], (μ_r, μ_β, μ_W)[1:n_dims]), interp_prices)
        etpf = extrapolate(itp, Line())
        sigmaitp[ttt] = Interpolations.scale(etpf, (r0, β0, logW0)[1:n_dims])
        
    end

    return E_ret, aE_ret, fwds, carry_ret, aEitp, sigmaitp
    
end

function reduceStates(state::AbstractVector{Float64}, mu::AbstractVector{Int})

    if all(y->y>0, mu)
        state_redux = state
    elseif mu[2] > 0
        state_redux = state[1:2]
    elseif mu[3] > 0
        state_redux = state[[1,3]]
    else
        state_redux = [state[1]]
    end

end

function reduceStates(state::AbstractVector{Int64}, mu::AbstractVector{Int})

    if all(y->y>0, mu)
        state_redux = state
    elseif mu[2] > 0
        state_redux = state[1:2]
    elseif mu[3] > 0
        state_redux = state[[1,3]]
    else
        state_redux = [state[1]]
    end

end


