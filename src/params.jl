# params.jl 
# contains the structures and functions parameterizing the model
# also contains all structures used to store results 
# authors: Rohan Kekre, Moritz Lenel, Federico Mainardi - 2023
# see https://github.com/KekreLenel/term for updates

# struct to store parameters, default values of benchmark parameterization
@with_kw struct params{R<:Real,I<:Integer}
    
    # Model parameters 
    # ------------------------------------------------------------------

    σ_r::R              = 0.0068          # std interest rate shock
    r_bar::R            = -0.0004         # avr. interest rate # -0.006
    κ_r::R              = 0.03            # ar coeff r process 0.65
    W_bar::R            = 0.002           # arb wealth
    W_center::R         = 0.07            # Wealth grid center: try to set as low as possible to limit truncation
    γ::R                = 2.0             # risk aversion
    τ_max::R            = 30.             # max maturity
    α::R                = 5.0             # PH demand slope 5.21
    δ_al::R             = 1.0             # demand shock short mat.
    θ::R                = 1.0             # PH demand
    θ0::R               = 1.0
    δ_th::Vector{R}     = [1.0, 10000.0]  # demand shock loadings
    δ_th0::Vector{R}    = [1.0, 10000.0]  # base demand loadings
    ξ::R                = 0.05            # wealth reversal 0.01
    σ_β::R              = 0.45            # std demand shock
    κ_β::R              = .08
    β_bar::R            = 0
    κ_m::R             = 0.25         # ar coeff monetary shock
    endo_wealth::I      = 1

    σ2_r::R       = σ_r^2              # var interest rate shock
    σ_r_unc::R    = sqrt(σ2_r/(2*κ_r)) # unconditional std dev
    σ2_β::R       = σ_β^2              # var demand shock
    σ_β_unc::R    = sqrt(σ2_β/(2*κ_β)) # unconditional std dev
   
    # QE parameters
    gdp::R = 14.45*1e+3 # Billion $

    # r* experiment
    r_bar_start::R = 0.0263539639810759
    r_bar_end::R   = 0.00389516292256526
    
    # Numerical parameters 
    # ------------------------------------------------------------------

    convergeCrit = 1E-05 # convergence criterion
    
    interp_prices = BSpline(Cubic(Line(OnGrid()))) # interpolation method for bond prices 
    interp_drifts = BSpline(Cubic(Line(OnGrid()))) # interpolation method for drifts and loadings
    
    # time steps 
    f_τ::I = 12 # number of periods per year for stored prices
    n_t::I = 6 # number of time steps per period for simulation
    n_τ::I = Int(round(τ_max*f_τ))+1
    τ::LinRange{R, I} = LinRange(0., τ_max, n_τ)
    Δτ::R    = (τ[end]-τ[1])/(n_τ -1)
    Δt::R    = Δτ/n_t

    # annual time steps
    y_idxs::Vector{I} = Int.(round.((1:τ_max)./Δτ) .+ 1)
    τ_y::Vector{R}   = τ[y_idxs]
    n_idxs::I = size(y_idxs,1)
    Δτ_y::I  = Int(τ_y[2]-τ_y[1])
    Δt_y::I  = Int(Δτ_y/Δτ)
    
    # tensor grid
    n_dims::I = 3
    μ_r::I = 5    
    μ_β::I = 5
    μ_W::I = 5
    mu::Vector{I}     = [μ_r, μ_β, μ_W]
    n_grid::I        = μ_r*μ_β*μ_W
    r_grd_wdth_L::R = 3.0  # witdh of r grid: x*std_dev of r process
    r_grd_wdth_U::R = 5.0  # witdh of r grid: x*std_dev of r process
    b_grd_wdth_L::R = 3.0  # witdh of β grid: x*std_dev of β process
    b_grd_wdth_U::R = 3.0  # witdh of β grid: x*std_dev of β process
    W_grd_wdth_L::R = 1.5              
    W_grd_wdth_U::R = 1.5             
    rbnds::Vector{R} =  [ (r_bar-r_grd_wdth_L*σ_r_unc*(mu[1]>1)), (r_bar+r_grd_wdth_U*σ_r_unc*(mu[1]>1))     ] 
    βbnds::Vector{R} =  [ (β_bar-b_grd_wdth_L*σ_β_unc*(mu[2]>1)), (β_bar+b_grd_wdth_U*σ_β_unc*(mu[2]>1))     ] 
    Wbnds::Vector{R} =  [ (log(W_center) - W_grd_wdth_L*(mu[3]>1)),    (log(W_center)  + W_grd_wdth_U*(mu[3]>1)) ]
    extrap::R        = 1.3
    rbnds_x::Vector{R} =  [ (r_bar-extrap*r_grd_wdth_L*σ_r_unc*(mu[1]>1)), (r_bar+extrap*r_grd_wdth_U*σ_r_unc*(mu[1]>1))     ] 
    βbnds_x::Vector{R} =  [ (β_bar-extrap*b_grd_wdth_L*σ_β_unc*(mu[2]>1)), (β_bar+extrap*b_grd_wdth_U*σ_β_unc*(mu[2]>1))     ] 
    Wbnds_x::Vector{R} =  [ (log(W_center) - W_grd_wdth_L*(mu[3]>1)),    (log(W_center)  + W_grd_wdth_U*(mu[3]>1)) ] 

    n_q::I = 3  # quadrature nodes per dimension
    
    adj_coeff::R = 3.0              
    
    # monte carlo runs
    n_mc::I          = 500 # number of MC sims
    
    # simulation parameters
    T_burn::I  = Int(100.0/Δτ) # number of periods to burn
    NN_sim::I  = 10  # number of long simulations (needs to be divisible by 2)
    N_sim::I   = 300 # number of short simulations within each long simulation
    Len_sim::I = Int(13/Δτ)+1 # length of short simulation in number of periods
    T_sim::I   = N_sim*Len_sim + T_burn # total number of periods per long simulation
    T_ss::I    = Int(1000.0/Δτ)
    T_paste::I = Int(25.0/Δτ)

    # IR parameters
    N_sample::I = 100
    T_ir::I  = Int(60/Δτ) 
    T_irr::I = Int(12*13+1)
    t_max::I = Int(5/Δτ)

    # Create bond demand vectors
    demand_τmax = τ_max 
    θvec  = demandfun(τ, θ, δ_th, demand_τmax) 
    θ0vec = demandfun(τ, θ0, δ_th0, demand_τmax)
    αvec  = alphafun(τ, θ0vec, α, δ_al, Δτ)
   
    # r* path  
    r_bar_shift = LinRange(r_bar_start, r_bar_end, T_irr)[:,1]

    # Program details 
    # ------------------------------------------------------------------
    calccarry::Bool        = true  # calculate expected carry returns
    calc_IRM::Bool         = false # solve for monetary shock
    calc_QE::Bool          = false # solve for QE shock
    calccarry_shocks::Bool = false # calculate carry returns for monetary or QE shocks
    calc_aE_sim::Bool      = false # calculate expected returns during simulation

    shock_size = 0.2
    cutoff = 0.05
end

function getParams()
    
    runall = true 
    # benchmark parameterization 
    prmBM1 = params{Float64, Int64}(; calc_IRM=runall, calc_QE=runall)
    # exogenous wealth, wealth calibrated to match slope of the yield curve in benchmark
    prmBM2 = params{Float64, Int64}(; adj_coeff=1.0E8, θ = 1.0-1.0/(1.0 + 3.0), n_dims = 2, μ_W = 1, ξ = 5.0, endo_wealth = 0, γ = 0.116, W_center = 0.002, calc_IRM=runall, calc_QE=runall) 
    # no demand shocks
    prmBM3 = params{Float64, Int64}(; σ_β = 1E-06)
    # lower duration 
    prmBM4 = params{Float64, Int64}(; W_center = 0.105, W_bar = 0.04, γ= 4.4, calc_IRM=runall)
    # lower demand elasticity
    prmBM5 = params{Float64, Int64}(; α = 0.0, γ=1.6, W_center = 0.055,  calc_IRM=runall, calc_QE=runall)
    # exogenous wealth, wealth calibrated to match average wealth level in benchmark
    prmBM6 = params{Float64, Int64}(; adj_coeff=1.0E8, μ_W = 1, ξ = 5.0, n_dims = 2, endo_wealth = 0, W_bar = 0.052, W_center=0.052)
    # less persistent monetary shock 
    prmBM7 = params{Float64, Int64}(; calc_IRM=runall, κ_m=1.0)
    
    allprms = [ 
                prmBM1,        
                prmBM2,         
                prmBM3,        
                prmBM4,        
                prmBM5,        
                prmBM6,        
                prmBM7,        
              ]
    
    return allprms

end

function getDictionary()

    calibs = Dict("benchmark" => 1, "exogenous_wealth" => 2, "no_demand_shocks" => 3, 
    "low_duration" => 4, "zero_alpha" => 5, "benchmark_exogenous_wealth" => 6, "high_kappa" => 7)

    return calibs

end


@with_kw struct FigParams{R<:Real}
    
    xgridvisible::Bool = false
    ygridvisible::Bool = false
    framevisible::Bool = true

    # fontsizes
    bigfont::R      = 12.0
    smallfont::R    = 10.0
    titlesize::R      = bigfont 
    yticklabelsize::R = bigfont  
    xticklabelsize::R = bigfont 
    ylabelsize::R     = bigfont 
    xlabelsize::R     = bigfont 
    labelsize::R      = smallfont
    
    # linewidth
    boxwidth = 0.5
    spinewidth::R = boxwidth
    xtickwidth::R = boxwidth
    ytickwidth::R = boxwidth
    framewidth::R = boxwidth
    spindwidth::R = boxwidth
    linewidth::R = 1.5

    # legend position
    position=:rt
    
    # legend padding
    padding=(3.0f0, 3.0f0, 3.0f0, 3.0f0)
    pt_per_unit::R=1.0

    strokewidth::R = 0.0
    markersize::R = 8.0
   
    # figure formats 
    rescaling::R = 0.85

    # arg lists 
    axisargs = Dict(:xgridvisible => xgridvisible, :ygridvisible => ygridvisible, :titlesize => titlesize, :yticklabelsize => yticklabelsize, :titlefont=>:regular, 
                :xticklabelsize => xticklabelsize, :xlabelsize => xlabelsize, :ylabelsize => ylabelsize, :spinewidth=>spinewidth, :xtickwidth=>xtickwidth, :ytickwidth=>ytickwidth)
    legargs = Dict(:labelsize=>labelsize, :labelfont=>:regular, :framevisible=>true, :framewidth=>framewidth, :position=>position,  :padding=>padding)
    format1 = Dict(:width=>0.9*rescaling*5*72, :height=>0.9*rescaling*5*72/1.618)
    format2 = Dict(:width=>rescaling*2.5*72, :height=>rescaling*2.5*72)
    format3 = Dict(:width=>rescaling*6.3*72, :height=>rescaling*1.75*72)
    format4 = Dict(:width=>0.85*rescaling*1.75*72, :height=>0.85*rescaling*1.75*72)
end


function demandfun(x,θ,δ,tmax)

    y = zeros(size(x))
    y[ x .<= tmax ] .= θ*(exp.(-δ[1]*x[ x .<= tmax ]) - exp.(-δ[2]*x[ x .<= tmax ]))

    return y
        
end

function alphafun(x,θ0,α,δ_al,Δτ)

    y = zeros(size(x))
    y[ θ0 .>0 ] .= α*exp.(-δ_al*x[ θ0 .>0 ])
    return y
        
end

import Base: +, -, *, /, ≈

# Type for multi-valued interpolation
struct MyVec{T}
    myvec::SVector{6,T}
end

# Interface for multi-valued interpolation
(+)(p1::MyVec, p2::MyVec) = MyVec(p1.myvec .+ p2.myvec)
(-)(p1::MyVec, p2::MyVec) = MyVec(p1.myvec .- p2.myvec)
(*)(n::Number, p::MyVec)  = MyVec(n.*p.myvec)
(*)(p::MyVec, n::Number)  = n*p
(/)(p::MyVec, n::Number)  = MyVec(p.myvec./n)
Base.zero(::Type{MyVec{T}}) where {T} = MyVec(zeros(SVector{6}))
Base.promote_rule(::Type{MyVec{T1}}, ::Type{T2}) where {T1,T2<:Number} = MyVec{promote_type(T1,T2)}
≈(p1::MyVec, p2::MyVec) = (p1.myvec ≈ p2.myvec)
          
# Type for multi-valued interpolation
struct MyVec2{T}
    myvec::SVector{361,T}
end

# Interface for multi-valued interpolation
(+)(p1::MyVec2, p2::MyVec2) = MyVec2(p1.myvec .+ p2.myvec)
(-)(p1::MyVec2, p2::MyVec2) = MyVec2(p1.myvec .- p2.myvec)
(*)(n::Number, p::MyVec2)  = MyVec2(n.*p.myvec)
(*)(p::MyVec2, n::Number)  = n*p
(/)(p::MyVec2, n::Number)  = MyVec2(p.myvec./n)
Base.zero(::Type{MyVec2{T}}) where {T} = MyVec2(zeros(SVector{361}))
Base.promote_rule(::Type{MyVec2{T1}}, ::Type{T2}) where {T1,T2<:Number} = MyVec2{promote_type(T1,T2)}
≈(p1::MyVec2, p2::MyVec2) = (p1.myvec ≈ p2.myvec)
          

# Prepare Inputs for Empirics

struct Empirics
    series::Vector{Float64}
    series_lci::Vector{Float64}
    series_uci::Vector{Float64}
end

struct RunEmpirics
    FB_data::Empirics
    CS_data::Empirics
    df1_real::Empirics
    df1::Empirics
    dy1::Empirics
    ymscatter::Vector{String}
    ymscatter_hf::Vector{String}
    dyscatter::Vector{Float64}
    dyscatter_hf::Vector{Float64}
    dfscatter20::Vector{Float64}
    avg_ret::Vector{Float64}
    avg_ret_hf::Vector{Float64}
    tips_issues::Matrix{Any}
    qe_purchases::Matrix{Float64}
    tp_5yf5y::Vector{Union{Missing, Float64}}
    log_dur_pd::Vector{Union{Missing, Float64}}
    aggincgap::Vector{Union{Missing, Float64}}
    yduration::Vector{String}
    qe1_10ye::Vector{Float64}
    yqqe::Vector{String}
end

# Prepare Inputs for Model

struct Stochsteady
    s::Vector{Float64} # steady state
    y::Vector{Float64} # yields
    f::Vector{Float64} # forwards
    E::Vector{Float64} # expected log returns
    ws::Float64
end

struct Impulseresponse
    r::Vector{Float64}
    β::Vector{Float64}
    W::Vector{Float64} 
    y::Array{Float64} 
    aE::Array{Float64}
    fwd::Array{Float64} 
    ret::Array{Float64}
    mret::Array{Float64}  
    CP::Vector{Array{Float64}}
    CPE::Vector{Array{Float64}}
    tp::Array{Float64}
    tp5f5::Vector{Float64}
    x::Array{Float64}
end

struct StateDepReg
    tp5f5::Float64
    dur::Float64
    ig::Float64
    tp5f5_restr::Float64
    ig_restr::Float64
end

struct Simulation
    mom::Vector{Float64} # moments
    FB::Vector{Float64}   # FB reg coeffs
    CS::Vector{Float64}   # CS reg coeffs
    y_avrg::Vector{Float64} # average yields
    fwd_avrg::Vector{Float64} # average forwards
    r::Vector{Float64}
    β::Vector{Float64}
    W::Vector{Float64}
    CP::Array{Float64}
    CPE::Array{Float64}
end

struct SimulationQ
    r::Vector{Float64}
    β::Vector{Float64}
    W::Vector{Float64}
end

struct Run
    ss::Stochsteady
    irM_sdreg::StateDepReg
    sim::Simulation
    simQ::SimulationQ
    irr::Impulseresponse
    irβ::Impulseresponse
    irW::Impulseresponse
    irM::Impulseresponse
    irQE::Impulseresponse
    irQE_lowW::Impulseresponse
    irrstar::Impulseresponse
    prm::params
end

