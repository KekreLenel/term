function makeCalibTable(moments::Vector{R}, FB::Vector{R}, ir::Impulseresponse, prms::params{R,I}, filename::String) where {R <: Real, I <: Integer}


    @unpack κ_r, σ_r, ξ, κ_β, σ_β, α, δ_th, r_bar, W_bar, γ = prms

    df10 = normalize1(vec(ir.fwd[9,:]), 2)
    df20 = normalize1(vec(ir.fwd[19,:]), 2)
    targets = [0.06, 1.54, 0.74, 0.56, 0.68, 20, 0.5, -0.46, -0.25]
    open(filename, "w") do f
        @printf(f, "\\begin{table}[htbp]  \n")
        @printf(f, "\\centering  \n")
        @printf(f, "\\bgroup  \n")
        @printf(f, "\\def\\arraystretch{1.25}  \n")
        @printf(f, "\\begin{tabular}{clclcc} \\hline  \n")
        @printf(f, "& Description & Value & Moment & Target & Model \\\\ \\hline  \n")
        @printf(f, "\\multicolumn{6}{l}{\\emph{Unconditional moments of yields and volumes}} \\\\  \n")
        @printf(f, "\$\\bar{r}\$  & mean short rate & %.4f  & \$y_{t}^{(1)}\$ & %.2f\\%% & %.2f\\%%   \\\\ \n", r_bar, targets[1], moments[1])
        @printf(f, "\$\\gamma\$  & arb. risk aversion & %.2f  & \$y_{t}^{(20)}-y_{t}^{(1)}\$ & %.2f\\%% & %.2f\\%%   \\\\ \n", γ, targets[2], moments[41])
        @printf(f, "\$\\sigma_r\$ & std. dev. short rate & %.3f & \$\\sigma ( y_{t}^{(20)} )\$  & %.2f\\%% & %.2f\\%% \\\\ \n", σ_r, targets[3], moments[27])
        @printf(f, "\$\\kappa_r\$ & mean rev. short rate & %.3f & \$\\sigma (\\Delta y_{t}^{(20)})\$ & %.2f\\%% & %.2f\\%% \\\\ \n", κ_r, targets[4], moments[28])
        @printf(f, "\$\\sigma_{\\beta}\$ & std. dev. demand & %.2f & \$ \\beta^{(10)}_{FB} \$ & %.2f & %.2f \\\\ \n", σ_β, targets[5], FB[9])
        @printf(f, "\$\\bar{W}\$ & arb.endowment & %.3f & duration & %.2f & %.2f \\\\ \n", W_bar, targets[6], moments[88])
        @printf(f, "\$\\kappa_{\\beta}\$ & mean rev. demand & %.2f & \$\\sigma (\\log \$ (wealth duration) ) & %.2f & %.2f \\\\ \n", κ_β, targets[7], moments[102])
        @printf(f, "\$\\alpha\$ & level price elast. & %.2f & \$ df_t^{(9,10)} \$ & %.2f\\%% & %.2f\\%% \\\\ \n", α, targets[8], only(df10*100))
        @printf(f, "\$\\xi\$ & persistence arb. wealth & %.2f & \$ df_t^{(19,20)} \$ & %.2f & %.2f \\\\ \n", ξ, targets[9], only(df20*100))
        @printf(f, "\\end{tabular} \n")
        @printf(f, "\\egroup \n")
        @printf(f, "\\caption{baseline calibration} \n")
        @printf(f, "\\label{tab:cal} \n")
        @printf(f, "\\end{table} \n")
    end

nothing

end

function makeDurTable(irM_sdreg::StateDepReg, filename::String)

    open(filename, "w") do f

        @printf(f, "\\begin{table}[tb]\\centering \n")
        @printf(f, "\\begin{tabular}{lC{18mm}C{18mm}C{18mm}} \\hline \n")
        @printf(f, "& \\multicolumn{3}{c}{Proxy for arb duration} \\\\  \n")
        @printf(f, "& 5-yr fwd, 5-yr TP & Log dealer dur. & \$-\$Dealer income gap \\\\ \\hline  \n")
        @printf(f, "Data & [0.09,0.91] & [0.09,0.59] & [-0.8,4.8] \\\\  \n")
        @printf(f, "Model & %.2f & %.2f & %.2f  \\\\ \n", irM_sdreg.tp5f5, irM_sdreg.dur, irM_sdreg.ig)
        @printf(f, "\\end{tabular} \n")
        @printf(f, "\\caption{\$\\Delta f_t^{(20)}\$ on \$\\Delta y_t^{(1)}\$, duration of arbitrageurs, and interaction given monetary shock: model vs. data} \n")
        @printf(f, "\\floatfoot{Notes: empirical estimates correspond to 90\\%% confidence interval from baseline estimates in Table \\ref{tab:statedepmp}.} \n")
        @printf(f, "\\label{tab:statedepmp_modelvsdata} \n")
        @printf(f, "\\end{table} \n ")
    end

nothing

end

function makeLongYieldsTable(moments1::Vector{R}, moments2::Vector{R}, filename::String) where {R <: Real}

    open(filename, "w") do f

        @printf(f, "\\begin{table}[htbp]\\centering  \n \\bgroup  \n \\def\\arraystretch{1.25}")
        @printf(f, "\\begin{tabular}{lcc} \\hline \n")
        @printf(f, "Moment & Model & \$\\xi \\rightarrow \\infty\$ \\\\ \\hline \n")
        @printf(f, "\$\\sigma(y_t^{(20)})\$ & %.3f\\%% & %.3f\\%% \\\\ \n", moments1[27], moments2[27])
        @printf(f, "\$\\sigma(\\sigma_{t-1}(y_t^{(20)}))\$ & %.3f\\%% & %.3f\\%% \\\\ \n", moments1[107], moments2[107]) 
        @printf(f, "\$y_t^{(20)}-y_t^{(1)}\$ & %.3f\\%% & %.3f\\%% \\\\ \\hline \n", moments1[41], moments2[41])
        @printf(f, "\\end{tabular} \n \\egroup \n")
        @printf(f, "\\caption{unconditional moments of long yields} \n")
        @printf(f, "\\floatfoot{Notes: \$\\sigma\$ denotes monthly standard deviation and last row is simple time-series average.} \n")
        @printf(f, "\\end{table} \n")
    end 
        
end

function makeAnalyticalTable(moments::Vector{R}, FB::Vector{R}, prms::params{R,I}, filename::String) where {R <: Real, I <: Integer}

    @unpack n_idxs, r_bar, β_bar = prms

    amoments = zeros(10)
    aterm = zeros(n_idxs)
    try
        amoments, aterm = calcAnalytical(prms, r_bar, β_bar)
    catch
        println("Analytical solution has no equilibrium")
    end

    open(filename, "w") do f
        @printf(f, "\\begin{table}[htbp]\\centering  \n \\bgroup  \n \\def\\arraystretch{1.25}")
        @printf(f, "\\caption{Numerical vs. analytical} \n")
        @printf(f, "\\begin{tabular}{llrr} \n")
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\multicolumn{3}{c}{} \\\\ \n")
        @printf(f, "Moment & Numerical  & Closed form \\\\ \n")
        @printf(f, "\\hline \n")
        @printf(f, "\$y_{t}^{(1)}\$ & %.2f & %.2f  \\\\ \n", moments[116], amoments[1])
        @printf(f, "\$y_{t}^{(10)}-y_{t}^{(1)}\$ & %.2f &  %.2f  \\\\ \n", moments[117], amoments[2])
        @printf(f, "\$y_{t}^{(20)}-y_{t}^{(1)}\$ & %.2f &  %.2f  \\\\ \n", moments[118], amoments[3])
        @printf(f, "\$ \\beta^{(10)}_{FB}\$ & %.2f & %.2f \\\\ \n", FB[9], amoments[4])
        @printf(f, "\$\\sigma(y_{t}^{(5)})\$ & %.2f & %.2f  \\\\ \n", moments[119], amoments[5])
        @printf(f, "\$\\sigma(y_{t}^{(10)})\$ & %.2f & %.2f  \\\\ \n", moments[120], amoments[6])
        @printf(f, "\$\\sigma(y_{t}^{(20)})\$ & %.2f & %.2f  \\\\ \n", moments[121], amoments[7])
        @printf(f, "\$\\sigma(\\Delta y_{t}^{(20)})\$ & %.2f & %.2f  \\\\ \n", moments[122], amoments[8])
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\end{tabular} \n \\egroup \n")
        @printf(f, "\\end{table} \n")
    end 
        
end

# additional tables
# ====================================================================================================================================


function makeNumericalTable(prms::params{R,I}, filename::String) where {R <: Real, I <: Integer}

    @unpack_params prms

    open(filename,"w") do f
        @printf(f, "\\begin{table}[htbp]\\centering \n")
        @printf(f, "\\caption{Numerical Inputs} \n")
        @printf(f, "\\begin{tabular}{l*{1}{ccc}} \n")
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\multicolumn{3}{c}{} \\\\ \n")
        @printf(f, "Parameter & Description & Value \\\\ \n")
        @printf(f, "\\hline \n")
        @printf(f, "\$\\mu_{r}\$ & Number of Points on \$r\$ Grid & %6.0f \\\\ \n", μ_r)
        @printf(f, "\$\\mu_{\\beta}\$ & Number of Points on \$\\beta\$ Grid & %6.0f \\\\ \n", μ_β)
        @printf(f, "\$\\mu_{W}\$ & Number of Points on \$W\$ Grid & %6.0f \\\\ \n", μ_W)
        @printf(f, "\$r_{width}\$ & Width of \$r\$ Grid & \$[%6.2f \\ %6.2f]*std(r)\$ \\\\ \n", r_grd_wdth_L, r_grd_wdth_U)
        @printf(f, "\$\\beta_{width}\$ & Width of \$\\beta\$ Grid & \$[%6.2f \\ %6.2f]*std(\\beta)\$ \\\\ \n", b_grd_wdth_L, b_grd_wdth_U)
        @printf(f, "\$W_{width}\$ & Width of \$\\log W\$ Grid & \$[%6.2f \\ %6.2f] \\\\ \n", W_grd_wdth_L, W_grd_wdth_U)
        @printf(f, "\$W_{center}\$ & center of W grid & %6.2f  \\\\ \n", W_center)
        @printf(f, "\$T_{sim}\$ & Number of simulated periods & %6.0f \\\\ \n", T_sim)
        @printf(f, "\$T_{ir}\$ & Time horizon for IR & %6.0f \\\\ \n", t_max)
        @printf(f, "\$N_{\\tau}\$ & Number of points on maturity grid & %6.0f \\\\ \n", n_τ)
        @printf(f, "\$\\bar{\\tau}\$ & Maximum maturity & %6.0f \\\\ \n", τ_max)
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\end{tabular} \n")
        @printf(f, "\\end{table} \n")
    end

nothing

end

function makeAddlTable(results::Run,  prms::params{R,I}, filename) where {R <: Real, I <: Integer}

    @unpack κ_m, adj_coeff, W_center, y_idxs = prms

    df2_m  = normalize1(vec(results.irM.fwd[1,:]), 2)
    df10_m = normalize1(vec(results.irM.fwd[9,:]), 2)
    df20_m = normalize1(vec(results.irM.fwd[19,:]), 2)
    dy1_m  = normalize1(vec(results.irM.y[y_idxs[1],:]), 2)
    dy10_m = normalize1(vec(results.irM.y[y_idxs[10],:]), 2)
    dy20_m = normalize1(vec(results.irM.y[y_idxs[20],:]), 2)
    dW_m     = normalize1(log.(results.irM.W), 2)

    df2_QE  = normalize1(vec(results.irQE_lowW.fwd[1,:]), 2)
    df10_QE = normalize1(vec(results.irQE_lowW.fwd[9,:]), 2)
    df20_QE = normalize1(vec(results.irQE_lowW.fwd[19,:]), 2)
    dy1_QE  = normalize1(vec(results.irQE_lowW.y[y_idxs[1],:]), 2)
    dy10_QE = normalize1(vec(results.irQE_lowW.y[y_idxs[10],:]), 2)
    dy20_QE = normalize1(vec(results.irQE_lowW.y[y_idxs[20],:]), 2)
    dW_QE     = normalize1(log.(results.irQE_lowW.W), 2)

    df2_QEavrg  = normalize1(vec(results.irQE.fwd[1,:]), 2)
    df10_QEavrg = normalize1(vec(results.irQE.fwd[9,:]), 2)
    df20_QEavrg = normalize1(vec(results.irQE.fwd[19,:]), 2)
    dy1_QEavrg  = normalize1(vec(results.irQE.y[y_idxs[1],:]), 2)
    dy10_QEavrg = normalize1(vec(results.irQE.y[y_idxs[10],:]), 2)
    dy20_QEavrg = normalize1(vec(results.irQE.y[y_idxs[20],:]), 2)
    dW_QEavrg     = normalize1(log.(results.irQE.W), 2)
    
    moments = results.sim.mom

    open(filename, "w") do f
        @printf(f, "\\begin{table}[htbp]\\centering  \n \\bgroup  \n \\def\\arraystretch{1.25}")
        @printf(f, "\\caption{Additional moments and parameters} \n")
        @printf(f, "\\begin{tabular}{llr} \n")
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\multicolumn{3}{c}{} \\\\ \n")
        @printf(f, "Moment & Description & Value  \\\\ \n")
        @printf(f, "\$\\kappa_m\$ & persistence mon. shock  & %.2f   \\\\ \n", κ_m)
        @printf(f, "\$b*W_{center}\$ & \$\\beta\$ adjustment  & %.2f   \\\\ \n", adj_coeff)
        @printf(f, "\\hline \n")
        @printf(f, "\\multicolumn{3}{c}{Monetary shock} \\\\ \n")
        @printf(f, "\$df^{(2)}/dy^{(1)} \$ & effect on 2-year forward  & %.3f \\\\ \n",    only(df2_m/dy1_m))
        @printf(f, "\$df^{(10)}/dy^{(1)} \$ & effect on 10-year forward  & %.3f \\\\ \n",  only(df10_m/dy1_m))
        @printf(f, "\$df^{(20)}/dy^{(1)} \$ & effect on 20-year forward  & %.3f \\\\ \n",  only(df20_m/dy1_m))
        @printf(f, "\$dy^{(10)}/dy^{(1)} \$ & effect on 10-year yield  & %.3f \\\\ \n",  only(dy10_m/dy1_m))
        @printf(f, "\$dy^{(20)}/dy^{(1)} \$ & effect on 20-year yield  & %.3f \\\\ \n",  only(dy20_m/dy1_m))
        @printf(f, "\$d\\log W\$ & effect on wealth  & %.1fbp \\\\ \n",         100*100*only(dW_m))
        @printf(f, "\\hline \n")
        @printf(f, "\\multicolumn{3}{c}{QE} \\\\ \n")
        @printf(f, "\$df^{(2)}\$ & effect on 2-year forward  & %.1fbp \\\\ \n",    100*100*only(df2_QE))
        @printf(f, "\$df^{(10)}\$ & effect on 10-year forward  & %.1fbp \\\\ \n",  100*100*only(df10_QE))
        @printf(f, "\$df^{(20)}\$ & effect on 20-year forward  & %.1fbp \\\\ \n",  100*100*only(df20_QE))
        @printf(f, "\$dy^{(10)}\$ & effect on 10-year yield  & %.1fbp \\\\ \n",    100*100*only(dy10_QE))
        @printf(f, "\$dy^{(20)}\$ & effect on 20-year yield  & %.1fbp \\\\ \n",    100*100*only(dy20_QE))
        @printf(f, "\$d\\log W\$ & effect on wealth  & %.1fbp \\\\ \n",            100*100*only(dW_QE))
        @printf(f, "\\hline \n")
        @printf(f, "\\multicolumn{3}{c}{QE (avrg W)} \\\\ \n")
        @printf(f, "\$df^{(2)} \$ & effect on 2-year forward  & %.1fbp \\\\ \n",    100*100*only( df2_QEavrg))
        @printf(f, "\$df^{(10)} \$ & effect on 10-year forward  & %.1fbp \\\\ \n",  100*100*only(df10_QEavrg))
        @printf(f, "\$df^{(20)} \$ & effect on 20-year forward  & %.1fbp \\\\ \n",  100*100*only(df20_QEavrg))
        @printf(f, "\$dy^{(10)} \$ & effect on 10-year yield  & %.1fbp \\\\ \n",    100*100*only(dy10_QEavrg))
        @printf(f, "\$dy^{(20)} \$ & effect on 20-year yield  & %.1fbp \\\\ \n",    100*100*only(dy20_QEavrg))
        @printf(f, "\$d\\log W\$ & effect on wealth  & %.1fbp \\\\ \n",             100*100*only(dW_QEavrg))
        @printf(f, "\\hline \n")
        @printf(f, "ss W & stoch. steady W & %.3f   \\\\ \n", moments[81])
        @printf(f, "\$W\$ & avrg. W & %.3f   \\\\ \n", moments[76])
        @printf(f, "\$\\sqrt{Var[\\log(W)]}\$ &  & %.2f   \\\\ \n", moments[77])
        @printf(f, "\$\\max \\log W - \\log W_{center} \$ &  & %.2f   \\\\ \n", moments[78])
        @printf(f, "\$\\min \\log W - \\log W_{center} \$ &  & %.2f   \\\\ \n", moments[79])
        @printf(f, "\$\\rho(W)\$ &  & %.2f   \\\\ \n", moments[80])
        @printf(f, "\\hline \n")
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\end{tabular} \n \\egroup \n")
        @printf(f, "\\end{table} \n")
    end 
        
end


function makeYieldMomsTable(moments::Vector{R}, filename::String) where {R <: Real}

    open(filename, "w") do f
    @printf(f, "\\begin{tabular}{|l*{5}{|cc|}} \n")
    @printf(f, "\\hline\\hline \n")
    @printf(f, "            & \\multicolumn{2}{|c|}{mean} & \\multicolumn{2}{|c|}{sigma} & \\multicolumn{2}{|c|}{sigma(d)} & \\multicolumn{2}{|c|}{rho} & \\multicolumn{2}{|c|}{rho(y1)} \\\\ \n")
    @printf(f, "\\hline \n")
    @printf(f, "            & data & model & data & model & data & model & data & model &  data & model \\\\ \n")
    @printf(f, "\\hline \n")
    @printf(f, "y1          & 0.06 & %.2f & 1.66 & %.2f & 1.75 & %.2f & 0.92 & %.2f &  1.00 & %.2f \\\\ \n", moments[1:5]... )
    @printf(f, "y2          & 0.10 & %.2f & 1.43 & %.2f & 1.28 & %.2f & 0.95 & %.2f &  0.97 & %.2f \\\\ \n", moments[6:10]... )
    @printf(f, "y3          & 0.24 & %.2f & 1.32 & %.2f & 1.07 & %.2f & 0.96 & %.2f &  0.93 & %.2f \\\\ \n", moments[11:15]... )
    @printf(f, "y5          & 0.55 & %.2f & 1.19 & %.2f & 0.88 & %.2f & 0.97 & %.2f &  0.86 & %.2f \\\\ \n", moments[16:20]... )
    @printf(f, "y10         & 1.14 & %.2f & 0.96 & %.2f & 0.69 & %.2f & 0.97 & %.2f &  0.75 & %.2f \\\\ \n", moments[21:25]... )
    @printf(f, "y20         & 1.60 & %.2f & 0.74 & %.2f & 0.56 & %.2f & 0.97 & %.2f &  0.61 & %.2f \\\\ \n", moments[26:30]... )
    @printf(f, "y20-y1      & 1.54 & %.2f & 1.34 & %.2f & 1.67 & %.2f & 0.89 & %.2f & -0.90 & %.2f \\\\ \n", moments[41:45]... )
    @printf(f, "y20-y2      & 1.50 & %.2f & 1.02 & %.2f & 1.18 & %.2f & 0.92 & %.2f & -0.91 & %.2f \\\\ \n", moments[46:50]... )
    @printf(f, "f1\\_2      & 0.15 & %.2f & 1.33 & %.2f & 0.97 & %.2f & 0.95 & %.2f &  0.84 & %.2f \\\\ \n", moments[51:55]... )
    @printf(f, "f1\\_5      & 1.15 & %.2f & 1.03 & %.2f & 0.80 & %.2f & 0.97 & %.2f &  0.63 & %.2f \\\\ \n", moments[56:60]... )
    @printf(f, "f1\\_10     & 1.98 & %.2f & 0.73 & %.2f & 0.62 & %.2f & 0.95 & %.2f &  0.42 & %.2f \\\\ \n", moments[61:65]... )
    @printf(f, "f1\\_20     & 1.90 & %.2f & 0.55 & %.2f & 0.59 & %.2f & 0.87 & %.2f &  0.29 & %.2f \\\\ \n", moments[66:70]... )
    @printf(f, "5y5 tp      &      & %.2f &      & %.2f &      & %.2f &      & %.2f &       & %.2f \\\\ \n", moments[108], moments[99], moments[109:111]...)
    @printf(f, "10y tp      &      & %.2f &      & %.2f &      & %.2f &      & %.2f &       & %.2f \\\\ \n", moments[112], moments[100], moments[113:115]...)
    @printf(f, "\\hline\\hline \n")
    @printf(f, "\\end{tabular} \n")
    end

end




function makeCPtable(results::Run, calib::String)

    # CP decomposition based on realized returns

    factor = 1/abs(results.irr.y[13,2] - results.irr.y[13,1])

    makeCPtable(results.irr.CP[2], results.irr.CP[1], factor, calib, "r","real")

    makeCPtable(results.irW.CP[2], results.irW.CP[1], factor, calib, "W","real")

    makeCPtable(results.irβ.CP[2], results.irβ.CP[1], factor, calib, "B","real")

    makeCPtable(results.irM.CP[2], results.irM.CP[1], factor, calib, "M","real")

    # CP decomposition based on expected returns

    factor = 1/abs(results.irr.y[13,2]-results.irr.y[13,1])
    
    makeCPtable(results.irr.CPE[2], results.irr.CPE[1], factor, calib, "r","exp")

    makeCPtable(results.irW.CPE[2], results.irW.CPE[1], factor, calib, "W","exp")

    makeCPtable(results.irβ.CPE[2], results.irβ.CPE[1], factor, calib, "B","exp")

    makeCPtable(results.irM.CPE[2], results.irM.CPE[1], factor, calib, "M","exp")


nothing

end


function makeCPtable(CPMatrix::Matrix{R}, CPMatrix_ss::Matrix{R}, factor::R, calib::String, state::String, method::String) where {R <: Real}

    CPMatrix = hcat(CPMatrix, zeros(size(CPMatrix,1),20))
    CPMatrix = vcat(CPMatrix, zeros(20,size(CPMatrix,2)))
    CPMatrix_ss = hcat(CPMatrix_ss, zeros(size(CPMatrix_ss,1),20))
    CPMatrix_ss = vcat(CPMatrix_ss, zeros(20,size(CPMatrix_ss,2)))

    pos = 1:1:19 
    CPMatrix = CPMatrix[pos,:]*100
    CPMatrix_ss = CPMatrix_ss[pos,:]*100

    if state == "r"
        dCPMatrix = (CPMatrix - CPMatrix_ss).*factor
    else
        dCPMatrix = (CPMatrix - CPMatrix_ss)*100
    end

    if state == "W"
        shock = "0.1 of SS"
    elseif state == "M"
        shock = "10 bps"
    else
        shock = "1 S.D."
    end

    open("../output/tables/CPDec_"*state*"_"*method*"_"*calib*".tex", "w") do f
        @printf(f, "\\begin{R}  \n")
        @printf(f, "\\begin{table}[htbp]\\centering  \n \\bgroup  \n \\def\\arraystretch{1.25}")
        @printf(f, "\\caption{Cochrane-Piazzesi Decomposition - %s Shock = %s} \n", state, shock)
        @printf(f, "\\tabcolsep=0.09cm\\begin{tabular}{l*{1}{cccccccccccccccccccccc}} \n")
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\multicolumn{22}{c}{} \\\\ \n")
        @printf(f, "\$\\tau\$ & Value & Sum & \$1\$ & \$2\$ & \$3\$ & \$4\$ & \$5\$ & \$6\$ & \$7\$ & \$8\$ & \$9\$ & \$10\$ & \$11\$ & \$12\$ & \$13\$ & \$14\$ & \$15\$ & \$16\$ & \$17\$ & \$18\$ & \$19\$  \\\\ \n")
        @printf(f, "\\hline \n")
        @printf(f, "\$2\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[1,1], sum(CPMatrix[1,2:end]), CPMatrix[1,2], CPMatrix[1,3], CPMatrix[1,4], CPMatrix[1,5], CPMatrix[1,6], CPMatrix[1,7], CPMatrix[1,8], CPMatrix[1,9], CPMatrix[1,10], CPMatrix[1,11], CPMatrix[1,12], CPMatrix[1,13], CPMatrix[1,14], CPMatrix[1,15], CPMatrix[1,16], CPMatrix[1,17], CPMatrix[1,18], CPMatrix[1,19], CPMatrix[1,20])
        @printf(f, "\$3\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[2,1], sum(CPMatrix[2,2:end]), CPMatrix[2,2], CPMatrix[2,3], CPMatrix[2,4], CPMatrix[2,5], CPMatrix[2,6], CPMatrix[2,7], CPMatrix[2,8], CPMatrix[2,9], CPMatrix[2,10], CPMatrix[2,11], CPMatrix[2,12], CPMatrix[2,13], CPMatrix[2,14], CPMatrix[2,15], CPMatrix[2,16], CPMatrix[2,17], CPMatrix[2,18], CPMatrix[2,19], CPMatrix[2,20])
        @printf(f, "\$4\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[3,1], sum(CPMatrix[3,2:end]), CPMatrix[3,2], CPMatrix[3,3], CPMatrix[3,4], CPMatrix[3,5], CPMatrix[3,6], CPMatrix[3,7], CPMatrix[3,8], CPMatrix[3,9], CPMatrix[3,10], CPMatrix[3,11], CPMatrix[3,12], CPMatrix[3,13], CPMatrix[3,14], CPMatrix[3,15], CPMatrix[3,16], CPMatrix[3,17], CPMatrix[3,18], CPMatrix[3,19], CPMatrix[3,20])
        @printf(f, "\$5\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[4,1], sum(CPMatrix[4,2:end]), CPMatrix[4,2], CPMatrix[4,3], CPMatrix[4,4], CPMatrix[4,5], CPMatrix[4,6], CPMatrix[4,7], CPMatrix[4,8], CPMatrix[4,9], CPMatrix[4,10], CPMatrix[4,11], CPMatrix[4,12], CPMatrix[4,13], CPMatrix[4,14], CPMatrix[4,15], CPMatrix[4,16], CPMatrix[4,17], CPMatrix[4,18], CPMatrix[4,19], CPMatrix[4,20])
        @printf(f, "\$6\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[5,1], sum(CPMatrix[5,2:end]), CPMatrix[5,2], CPMatrix[5,3], CPMatrix[5,4], CPMatrix[5,5], CPMatrix[5,6], CPMatrix[5,7], CPMatrix[5,8], CPMatrix[5,9], CPMatrix[5,10], CPMatrix[5,11], CPMatrix[5,12], CPMatrix[5,13], CPMatrix[5,14], CPMatrix[5,15], CPMatrix[5,16], CPMatrix[5,17], CPMatrix[5,18], CPMatrix[5,19], CPMatrix[5,20])
        @printf(f, "\$7\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[6,1], sum(CPMatrix[6,2:end]), CPMatrix[6,2], CPMatrix[6,3], CPMatrix[6,4], CPMatrix[6,5], CPMatrix[6,6], CPMatrix[6,7], CPMatrix[6,8], CPMatrix[6,9], CPMatrix[6,10], CPMatrix[6,11], CPMatrix[6,12], CPMatrix[6,13], CPMatrix[6,14], CPMatrix[6,15], CPMatrix[6,16], CPMatrix[6,17], CPMatrix[6,18], CPMatrix[6,19], CPMatrix[6,20])
        @printf(f, "\$8\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[7,1], sum(CPMatrix[7,2:end]), CPMatrix[7,2], CPMatrix[7,3], CPMatrix[7,4], CPMatrix[7,5], CPMatrix[7,6], CPMatrix[7,7], CPMatrix[7,8], CPMatrix[7,9], CPMatrix[7,10], CPMatrix[7,11], CPMatrix[7,12], CPMatrix[7,13], CPMatrix[7,14], CPMatrix[7,15], CPMatrix[7,16], CPMatrix[7,17], CPMatrix[7,18], CPMatrix[7,19], CPMatrix[7,20])
        @printf(f, "\$9\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[8,1], sum(CPMatrix[8,2:end]), CPMatrix[8,2], CPMatrix[8,3], CPMatrix[8,4], CPMatrix[8,5], CPMatrix[8,6], CPMatrix[8,7], CPMatrix[8,8], CPMatrix[8,9], CPMatrix[8,10], CPMatrix[8,11], CPMatrix[8,12], CPMatrix[8,13], CPMatrix[8,14], CPMatrix[8,15], CPMatrix[8,16], CPMatrix[8,17], CPMatrix[8,18], CPMatrix[8,19], CPMatrix[8,20])
        @printf(f, "\$10\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[9,1], sum(CPMatrix[9,2:end]), CPMatrix[9,2], CPMatrix[9,3], CPMatrix[9,4], CPMatrix[9,5], CPMatrix[9,6], CPMatrix[9,7], CPMatrix[9,8], CPMatrix[9,9], CPMatrix[9,10], CPMatrix[9,11], CPMatrix[9,12], CPMatrix[9,13], CPMatrix[9,14], CPMatrix[9,15], CPMatrix[9,16], CPMatrix[9,17], CPMatrix[9,18], CPMatrix[9,19], CPMatrix[9,20])
        @printf(f, "\$11\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[10,1], sum(CPMatrix[10,2:end]), CPMatrix[10,2], CPMatrix[10,3], CPMatrix[10,4], CPMatrix[10,5], CPMatrix[10,6], CPMatrix[10,7], CPMatrix[10,8], CPMatrix[10,9], CPMatrix[10,10], CPMatrix[10,11], CPMatrix[10,12], CPMatrix[10,13], CPMatrix[10,14], CPMatrix[10,15], CPMatrix[10,16], CPMatrix[10,17], CPMatrix[10,18], CPMatrix[10,19], CPMatrix[10,20])
        @printf(f, "\$12\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[11,1], sum(CPMatrix[11,2:end]), CPMatrix[11,2], CPMatrix[11,3], CPMatrix[11,4], CPMatrix[11,5], CPMatrix[11,6], CPMatrix[11,7], CPMatrix[11,8], CPMatrix[11,9], CPMatrix[11,10], CPMatrix[11,11], CPMatrix[11,12], CPMatrix[11,13], CPMatrix[11,14], CPMatrix[11,15], CPMatrix[11,16], CPMatrix[11,17], CPMatrix[11,18], CPMatrix[11,19], CPMatrix[11,20])
        @printf(f, "\$13\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[12,1], sum(CPMatrix[12,2:end]), CPMatrix[12,2], CPMatrix[12,3], CPMatrix[12,4], CPMatrix[12,5], CPMatrix[12,6], CPMatrix[12,7], CPMatrix[12,8], CPMatrix[12,9], CPMatrix[12,10], CPMatrix[12,11], CPMatrix[12,12], CPMatrix[12,13], CPMatrix[12,14], CPMatrix[12,15], CPMatrix[12,16], CPMatrix[12,17], CPMatrix[12,18], CPMatrix[12,19], CPMatrix[12,20])
        @printf(f, "\$14\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[13,1], sum(CPMatrix[13,2:end]), CPMatrix[13,2], CPMatrix[13,3], CPMatrix[13,4], CPMatrix[13,5], CPMatrix[13,6], CPMatrix[13,7], CPMatrix[13,8], CPMatrix[13,9], CPMatrix[13,10], CPMatrix[13,11], CPMatrix[13,12], CPMatrix[13,13], CPMatrix[13,14], CPMatrix[13,15], CPMatrix[13,16], CPMatrix[13,17], CPMatrix[13,18], CPMatrix[13,19], CPMatrix[13,20])
        @printf(f, "\$15\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[14,1], sum(CPMatrix[14,2:end]), CPMatrix[14,2], CPMatrix[14,3], CPMatrix[14,4], CPMatrix[14,5], CPMatrix[14,6], CPMatrix[14,7], CPMatrix[14,8], CPMatrix[14,9], CPMatrix[14,10], CPMatrix[14,11], CPMatrix[14,12], CPMatrix[14,13], CPMatrix[14,14], CPMatrix[14,15], CPMatrix[14,16], CPMatrix[14,17], CPMatrix[14,18], CPMatrix[14,19], CPMatrix[14,20])
        @printf(f, "\$16\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[15,1], sum(CPMatrix[15,2:end]), CPMatrix[15,2], CPMatrix[15,3], CPMatrix[15,4], CPMatrix[15,5], CPMatrix[15,6], CPMatrix[15,7], CPMatrix[15,8], CPMatrix[15,9], CPMatrix[15,10], CPMatrix[15,11], CPMatrix[15,12], CPMatrix[15,13], CPMatrix[15,14], CPMatrix[15,15], CPMatrix[15,16], CPMatrix[15,17], CPMatrix[15,18], CPMatrix[15,19], CPMatrix[15,20])
        @printf(f, "\$17\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[16,1], sum(CPMatrix[16,2:end]), CPMatrix[16,2], CPMatrix[16,3], CPMatrix[16,4], CPMatrix[16,5], CPMatrix[16,6], CPMatrix[16,7], CPMatrix[16,8], CPMatrix[16,9], CPMatrix[16,10], CPMatrix[16,11], CPMatrix[16,12], CPMatrix[16,13], CPMatrix[16,14], CPMatrix[16,15], CPMatrix[16,16], CPMatrix[16,17], CPMatrix[16,18], CPMatrix[16,19], CPMatrix[16,20])
        @printf(f, "\$18\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[17,1], sum(CPMatrix[17,2:end]), CPMatrix[17,2], CPMatrix[17,3], CPMatrix[17,4], CPMatrix[17,5], CPMatrix[17,6], CPMatrix[17,7], CPMatrix[17,8], CPMatrix[17,9], CPMatrix[17,10], CPMatrix[17,11], CPMatrix[17,12], CPMatrix[17,13], CPMatrix[17,14], CPMatrix[17,15], CPMatrix[17,16], CPMatrix[17,17], CPMatrix[17,18], CPMatrix[17,19], CPMatrix[17,20])
        @printf(f, "\$19\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[18,1], sum(CPMatrix[18,2:end]), CPMatrix[18,2], CPMatrix[18,3], CPMatrix[18,4], CPMatrix[18,5], CPMatrix[18,6], CPMatrix[18,7], CPMatrix[18,8], CPMatrix[18,9], CPMatrix[18,10], CPMatrix[18,11], CPMatrix[18,12], CPMatrix[18,13], CPMatrix[18,14], CPMatrix[18,15], CPMatrix[18,16], CPMatrix[18,17], CPMatrix[18,18], CPMatrix[18,19], CPMatrix[18,20])
        @printf(f, "\$20\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[19,1], sum(CPMatrix[19,2:end]), CPMatrix[19,2], CPMatrix[19,3], CPMatrix[19,4], CPMatrix[19,5], CPMatrix[19,6], CPMatrix[19,7], CPMatrix[19,8], CPMatrix[19,9], CPMatrix[19,10], CPMatrix[19,11], CPMatrix[19,12], CPMatrix[19,13], CPMatrix[19,14], CPMatrix[19,15], CPMatrix[19,16], CPMatrix[19,17], CPMatrix[19,18], CPMatrix[19,19], CPMatrix[19,20])
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\end{tabular} \n \\egroup \n")
        @printf(f, "\\end{table} \n")
        @printf(f, "\\end{landscape} \n")
    end 

    CPMatrix = CPMatrix_ss

    open("../output/tables/CPDec_ss_"*state*"_"*method*"_"*calib*".tex", "w") do f
        @printf(f, "\\begin{landscape}  \n")
        @printf(f, "\\begin{table}[htbp]\\centering  \n \\bgroup  \n \\def\\arraystretch{1.25}")
        @printf(f, "\\caption{Cochrane-Piazzesi Decomposition - SS} \n")
        @printf(f, "\\tabcolsep=0.09cm\\begin{tabular}{l*{1}{cccccccccccccccccccccc}} \n")
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\multicolumn{22}{c}{} \\\\ \n")
        @printf(f, "\$\\tau\$ & Value & Sum & \$1\$ & \$2\$ & \$3\$ & \$4\$ & \$5\$ & \$6\$ & \$7\$ & \$8\$ & \$9\$ & \$10\$ & \$11\$ & \$12\$ & \$13\$ & \$14\$ & \$15\$ & \$16\$ & \$17\$ & \$18\$ & \$19\$  \\\\ \n")
        @printf(f, "\\hline \n")
        @printf(f, "\$2\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[1,1], sum(CPMatrix[1,2:end]), CPMatrix[1,2], CPMatrix[1,3], CPMatrix[1,4], CPMatrix[1,5], CPMatrix[1,6], CPMatrix[1,7], CPMatrix[1,8], CPMatrix[1,9], CPMatrix[1,10], CPMatrix[1,11], CPMatrix[1,12], CPMatrix[1,13], CPMatrix[1,14], CPMatrix[1,15], CPMatrix[1,16], CPMatrix[1,17], CPMatrix[1,18], CPMatrix[1,19], CPMatrix[1,20])
        @printf(f, "\$3\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[2,1], sum(CPMatrix[2,2:end]), CPMatrix[2,2], CPMatrix[2,3], CPMatrix[2,4], CPMatrix[2,5], CPMatrix[2,6], CPMatrix[2,7], CPMatrix[2,8], CPMatrix[2,9], CPMatrix[2,10], CPMatrix[2,11], CPMatrix[2,12], CPMatrix[2,13], CPMatrix[2,14], CPMatrix[2,15], CPMatrix[2,16], CPMatrix[2,17], CPMatrix[2,18], CPMatrix[2,19], CPMatrix[2,20])
        @printf(f, "\$4\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[3,1], sum(CPMatrix[3,2:end]), CPMatrix[3,2], CPMatrix[3,3], CPMatrix[3,4], CPMatrix[3,5], CPMatrix[3,6], CPMatrix[3,7], CPMatrix[3,8], CPMatrix[3,9], CPMatrix[3,10], CPMatrix[3,11], CPMatrix[3,12], CPMatrix[3,13], CPMatrix[3,14], CPMatrix[3,15], CPMatrix[3,16], CPMatrix[3,17], CPMatrix[3,18], CPMatrix[3,19], CPMatrix[3,20])
        @printf(f, "\$5\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[4,1], sum(CPMatrix[4,2:end]), CPMatrix[4,2], CPMatrix[4,3], CPMatrix[4,4], CPMatrix[4,5], CPMatrix[4,6], CPMatrix[4,7], CPMatrix[4,8], CPMatrix[4,9], CPMatrix[4,10], CPMatrix[4,11], CPMatrix[4,12], CPMatrix[4,13], CPMatrix[4,14], CPMatrix[4,15], CPMatrix[4,16], CPMatrix[4,17], CPMatrix[4,18], CPMatrix[4,19], CPMatrix[4,20])
        @printf(f, "\$6\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[5,1], sum(CPMatrix[5,2:end]), CPMatrix[5,2], CPMatrix[5,3], CPMatrix[5,4], CPMatrix[5,5], CPMatrix[5,6], CPMatrix[5,7], CPMatrix[5,8], CPMatrix[5,9], CPMatrix[5,10], CPMatrix[5,11], CPMatrix[5,12], CPMatrix[5,13], CPMatrix[5,14], CPMatrix[5,15], CPMatrix[5,16], CPMatrix[5,17], CPMatrix[5,18], CPMatrix[5,19], CPMatrix[5,20])
        @printf(f, "\$7\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[6,1], sum(CPMatrix[6,2:end]), CPMatrix[6,2], CPMatrix[6,3], CPMatrix[6,4], CPMatrix[6,5], CPMatrix[6,6], CPMatrix[6,7], CPMatrix[6,8], CPMatrix[6,9], CPMatrix[6,10], CPMatrix[6,11], CPMatrix[6,12], CPMatrix[6,13], CPMatrix[6,14], CPMatrix[6,15], CPMatrix[6,16], CPMatrix[6,17], CPMatrix[6,18], CPMatrix[6,19], CPMatrix[6,20])
        @printf(f, "\$8\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[7,1], sum(CPMatrix[7,2:end]), CPMatrix[7,2], CPMatrix[7,3], CPMatrix[7,4], CPMatrix[7,5], CPMatrix[7,6], CPMatrix[7,7], CPMatrix[7,8], CPMatrix[7,9], CPMatrix[7,10], CPMatrix[7,11], CPMatrix[7,12], CPMatrix[7,13], CPMatrix[7,14], CPMatrix[7,15], CPMatrix[7,16], CPMatrix[7,17], CPMatrix[7,18], CPMatrix[7,19], CPMatrix[7,20])
        @printf(f, "\$9\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[8,1], sum(CPMatrix[8,2:end]), CPMatrix[8,2], CPMatrix[8,3], CPMatrix[8,4], CPMatrix[8,5], CPMatrix[8,6], CPMatrix[8,7], CPMatrix[8,8], CPMatrix[8,9], CPMatrix[8,10], CPMatrix[8,11], CPMatrix[8,12], CPMatrix[8,13], CPMatrix[8,14], CPMatrix[8,15], CPMatrix[8,16], CPMatrix[8,17], CPMatrix[8,18], CPMatrix[8,19], CPMatrix[8,20])
        @printf(f, "\$10\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[9,1], sum(CPMatrix[9,2:end]), CPMatrix[9,2], CPMatrix[9,3], CPMatrix[9,4], CPMatrix[9,5], CPMatrix[9,6], CPMatrix[9,7], CPMatrix[9,8], CPMatrix[9,9], CPMatrix[9,10], CPMatrix[9,11], CPMatrix[9,12], CPMatrix[9,13], CPMatrix[9,14], CPMatrix[9,15], CPMatrix[9,16], CPMatrix[9,17], CPMatrix[9,18], CPMatrix[9,19], CPMatrix[9,20])
        @printf(f, "\$11\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[10,1], sum(CPMatrix[10,2:end]), CPMatrix[10,2], CPMatrix[10,3], CPMatrix[10,4], CPMatrix[10,5], CPMatrix[10,6], CPMatrix[10,7], CPMatrix[10,8], CPMatrix[10,9], CPMatrix[10,10], CPMatrix[10,11], CPMatrix[10,12], CPMatrix[10,13], CPMatrix[10,14], CPMatrix[10,15], CPMatrix[10,16], CPMatrix[10,17], CPMatrix[10,18], CPMatrix[10,19], CPMatrix[10,20])
        @printf(f, "\$12\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[11,1], sum(CPMatrix[11,2:end]), CPMatrix[11,2], CPMatrix[11,3], CPMatrix[11,4], CPMatrix[11,5], CPMatrix[11,6], CPMatrix[11,7], CPMatrix[11,8], CPMatrix[11,9], CPMatrix[11,10], CPMatrix[11,11], CPMatrix[11,12], CPMatrix[11,13], CPMatrix[11,14], CPMatrix[11,15], CPMatrix[11,16], CPMatrix[11,17], CPMatrix[11,18], CPMatrix[11,19], CPMatrix[11,20])
        @printf(f, "\$13\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[12,1], sum(CPMatrix[12,2:end]), CPMatrix[12,2], CPMatrix[12,3], CPMatrix[12,4], CPMatrix[12,5], CPMatrix[12,6], CPMatrix[12,7], CPMatrix[12,8], CPMatrix[12,9], CPMatrix[12,10], CPMatrix[12,11], CPMatrix[12,12], CPMatrix[12,13], CPMatrix[12,14], CPMatrix[12,15], CPMatrix[12,16], CPMatrix[12,17], CPMatrix[12,18], CPMatrix[12,19], CPMatrix[12,20])
        @printf(f, "\$14\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[13,1], sum(CPMatrix[13,2:end]), CPMatrix[13,2], CPMatrix[13,3], CPMatrix[13,4], CPMatrix[13,5], CPMatrix[13,6], CPMatrix[13,7], CPMatrix[13,8], CPMatrix[13,9], CPMatrix[13,10], CPMatrix[13,11], CPMatrix[13,12], CPMatrix[13,13], CPMatrix[13,14], CPMatrix[13,15], CPMatrix[13,16], CPMatrix[13,17], CPMatrix[13,18], CPMatrix[13,19], CPMatrix[13,20])
        @printf(f, "\$15\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[14,1], sum(CPMatrix[14,2:end]), CPMatrix[14,2], CPMatrix[14,3], CPMatrix[14,4], CPMatrix[14,5], CPMatrix[14,6], CPMatrix[14,7], CPMatrix[14,8], CPMatrix[14,9], CPMatrix[14,10], CPMatrix[14,11], CPMatrix[14,12], CPMatrix[14,13], CPMatrix[14,14], CPMatrix[14,15], CPMatrix[14,16], CPMatrix[14,17], CPMatrix[14,18], CPMatrix[14,19], CPMatrix[14,20])
        @printf(f, "\$16\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[15,1], sum(CPMatrix[15,2:end]), CPMatrix[15,2], CPMatrix[15,3], CPMatrix[15,4], CPMatrix[15,5], CPMatrix[15,6], CPMatrix[15,7], CPMatrix[15,8], CPMatrix[15,9], CPMatrix[15,10], CPMatrix[15,11], CPMatrix[15,12], CPMatrix[15,13], CPMatrix[15,14], CPMatrix[15,15], CPMatrix[15,16], CPMatrix[15,17], CPMatrix[15,18], CPMatrix[15,19], CPMatrix[15,20])
        @printf(f, "\$17\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[16,1], sum(CPMatrix[16,2:end]), CPMatrix[16,2], CPMatrix[16,3], CPMatrix[16,4], CPMatrix[16,5], CPMatrix[16,6], CPMatrix[16,7], CPMatrix[16,8], CPMatrix[16,9], CPMatrix[16,10], CPMatrix[16,11], CPMatrix[16,12], CPMatrix[16,13], CPMatrix[16,14], CPMatrix[16,15], CPMatrix[16,16], CPMatrix[16,17], CPMatrix[16,18], CPMatrix[16,19], CPMatrix[16,20])
        @printf(f, "\$18\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[17,1], sum(CPMatrix[17,2:end]), CPMatrix[17,2], CPMatrix[17,3], CPMatrix[17,4], CPMatrix[17,5], CPMatrix[17,6], CPMatrix[17,7], CPMatrix[17,8], CPMatrix[17,9], CPMatrix[17,10], CPMatrix[17,11], CPMatrix[17,12], CPMatrix[17,13], CPMatrix[17,14], CPMatrix[17,15], CPMatrix[17,16], CPMatrix[17,17], CPMatrix[17,18], CPMatrix[17,19], CPMatrix[17,20])
        @printf(f, "\$19\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[18,1], sum(CPMatrix[18,2:end]), CPMatrix[18,2], CPMatrix[18,3], CPMatrix[18,4], CPMatrix[18,5], CPMatrix[18,6], CPMatrix[18,7], CPMatrix[18,8], CPMatrix[18,9], CPMatrix[18,10], CPMatrix[18,11], CPMatrix[18,12], CPMatrix[18,13], CPMatrix[18,14], CPMatrix[18,15], CPMatrix[18,16], CPMatrix[18,17], CPMatrix[18,18], CPMatrix[18,19], CPMatrix[18,20])
        @printf(f, "\$20\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[19,1], sum(CPMatrix[19,2:end]), CPMatrix[19,2], CPMatrix[19,3], CPMatrix[19,4], CPMatrix[19,5], CPMatrix[19,6], CPMatrix[19,7], CPMatrix[19,8], CPMatrix[19,9], CPMatrix[19,10], CPMatrix[19,11], CPMatrix[19,12], CPMatrix[19,13], CPMatrix[19,14], CPMatrix[19,15], CPMatrix[19,16], CPMatrix[19,17], CPMatrix[19,18], CPMatrix[19,19], CPMatrix[19,20])
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\end{tabular} \n \\egroup \n")
        @printf(f, "\\end{table} \n")
        @printf(f, "\\end{landscape} \n")
    end 

    CPMatrix = dCPMatrix

    open("../output/tables/CPDec_diff_"*state*"_"*method*"_"*calib*".tex", "w") do f
        @printf(f, "\\begin{landscape}  \n")
        @printf(f, "\\begin{table}[htbp]\\centering  \n \\bgroup  \n \\def\\arraystretch{1.25}")
        @printf(f, "\\caption{Cochrane-Piazzesi Decomposition - Difference %s} \n", state)
        @printf(f, "\\tabcolsep=0.09cm\\begin{tabular}{l*{1}{cccccccccccccccccccccc}} \n")
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\multicolumn{22}{c}{} \\\\ \n")
        @printf(f, "\$\\tau\$ & Value & Sum & \$1\$ & \$2\$ & \$3\$ & \$4\$ & \$5\$ & \$6\$ & \$7\$ & \$8\$ & \$9\$ & \$10\$ & \$11\$ & \$12\$ & \$13\$ & \$14\$ & \$15\$ & \$16\$ & \$17\$ & \$18\$ & \$19\$  \\\\ \n")
        @printf(f, "\\hline \n")
        @printf(f, "\$2\$  & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[1,1], sum(CPMatrix[1,2:end]), CPMatrix[1,2], CPMatrix[1,3], CPMatrix[1,4], CPMatrix[1,5], CPMatrix[1,6], CPMatrix[1,7], CPMatrix[1,8], CPMatrix[1,9], CPMatrix[1,10], CPMatrix[1,11], CPMatrix[1,12], CPMatrix[1,13], CPMatrix[1,14], CPMatrix[1,15], CPMatrix[1,16], CPMatrix[1,17], CPMatrix[1,18], CPMatrix[1,19], CPMatrix[1,20])
        @printf(f, "\$3\$  & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[2,1], sum(CPMatrix[2,2:end]), CPMatrix[2,2], CPMatrix[2,3], CPMatrix[2,4], CPMatrix[2,5], CPMatrix[2,6], CPMatrix[2,7], CPMatrix[2,8], CPMatrix[2,9], CPMatrix[2,10], CPMatrix[2,11], CPMatrix[2,12], CPMatrix[2,13], CPMatrix[2,14], CPMatrix[2,15], CPMatrix[2,16], CPMatrix[2,17], CPMatrix[2,18], CPMatrix[2,19], CPMatrix[2,20])
        @printf(f, "\$4\$  & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[3,1], sum(CPMatrix[3,2:end]), CPMatrix[3,2], CPMatrix[3,3], CPMatrix[3,4], CPMatrix[3,5], CPMatrix[3,6], CPMatrix[3,7], CPMatrix[3,8], CPMatrix[3,9], CPMatrix[3,10], CPMatrix[3,11], CPMatrix[3,12], CPMatrix[3,13], CPMatrix[3,14], CPMatrix[3,15], CPMatrix[3,16], CPMatrix[3,17], CPMatrix[3,18], CPMatrix[3,19], CPMatrix[3,20])
        @printf(f, "\$5\$  & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[4,1], sum(CPMatrix[4,2:end]), CPMatrix[4,2], CPMatrix[4,3], CPMatrix[4,4], CPMatrix[4,5], CPMatrix[4,6], CPMatrix[4,7], CPMatrix[4,8], CPMatrix[4,9], CPMatrix[4,10], CPMatrix[4,11], CPMatrix[4,12], CPMatrix[4,13], CPMatrix[4,14], CPMatrix[4,15], CPMatrix[4,16], CPMatrix[4,17], CPMatrix[4,18], CPMatrix[4,19], CPMatrix[4,20])
        @printf(f, "\$6\$  & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[5,1], sum(CPMatrix[5,2:end]), CPMatrix[5,2], CPMatrix[5,3], CPMatrix[5,4], CPMatrix[5,5], CPMatrix[5,6], CPMatrix[5,7], CPMatrix[5,8], CPMatrix[5,9], CPMatrix[5,10], CPMatrix[5,11], CPMatrix[5,12], CPMatrix[5,13], CPMatrix[5,14], CPMatrix[5,15], CPMatrix[5,16], CPMatrix[5,17], CPMatrix[5,18], CPMatrix[5,19], CPMatrix[5,20])
        @printf(f, "\$7\$  & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[6,1], sum(CPMatrix[6,2:end]), CPMatrix[6,2], CPMatrix[6,3], CPMatrix[6,4], CPMatrix[6,5], CPMatrix[6,6], CPMatrix[6,7], CPMatrix[6,8], CPMatrix[6,9], CPMatrix[6,10], CPMatrix[6,11], CPMatrix[6,12], CPMatrix[6,13], CPMatrix[6,14], CPMatrix[6,15], CPMatrix[6,16], CPMatrix[6,17], CPMatrix[6,18], CPMatrix[6,19], CPMatrix[6,20])
        @printf(f, "\$8\$  & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[7,1], sum(CPMatrix[7,2:end]), CPMatrix[7,2], CPMatrix[7,3], CPMatrix[7,4], CPMatrix[7,5], CPMatrix[7,6], CPMatrix[7,7], CPMatrix[7,8], CPMatrix[7,9], CPMatrix[7,10], CPMatrix[7,11], CPMatrix[7,12], CPMatrix[7,13], CPMatrix[7,14], CPMatrix[7,15], CPMatrix[7,16], CPMatrix[7,17], CPMatrix[7,18], CPMatrix[7,19], CPMatrix[7,20])
        @printf(f, "\$9\$  & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[8,1], sum(CPMatrix[8,2:end]), CPMatrix[8,2], CPMatrix[8,3], CPMatrix[8,4], CPMatrix[8,5], CPMatrix[8,6], CPMatrix[8,7], CPMatrix[8,8], CPMatrix[8,9], CPMatrix[8,10], CPMatrix[8,11], CPMatrix[8,12], CPMatrix[8,13], CPMatrix[8,14], CPMatrix[8,15], CPMatrix[8,16], CPMatrix[8,17], CPMatrix[8,18], CPMatrix[8,19], CPMatrix[8,20])
        @printf(f, "\$10\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[9,1], sum(CPMatrix[9,2:end]), CPMatrix[9,2], CPMatrix[9,3], CPMatrix[9,4], CPMatrix[9,5], CPMatrix[9,6], CPMatrix[9,7], CPMatrix[9,8], CPMatrix[9,9], CPMatrix[9,10], CPMatrix[9,11], CPMatrix[9,12], CPMatrix[9,13], CPMatrix[9,14], CPMatrix[9,15], CPMatrix[9,16], CPMatrix[9,17], CPMatrix[9,18], CPMatrix[9,19], CPMatrix[9,20])
        @printf(f, "\$11\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[10,1], sum(CPMatrix[10,2:end]), CPMatrix[10,2], CPMatrix[10,3], CPMatrix[10,4], CPMatrix[10,5], CPMatrix[10,6], CPMatrix[10,7], CPMatrix[10,8], CPMatrix[10,9], CPMatrix[10,10], CPMatrix[10,11], CPMatrix[10,12], CPMatrix[10,13], CPMatrix[10,14], CPMatrix[10,15], CPMatrix[10,16], CPMatrix[10,17], CPMatrix[10,18], CPMatrix[10,19], CPMatrix[10,20])
        @printf(f, "\$12\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[11,1], sum(CPMatrix[11,2:end]), CPMatrix[11,2], CPMatrix[11,3], CPMatrix[11,4], CPMatrix[11,5], CPMatrix[11,6], CPMatrix[11,7], CPMatrix[11,8], CPMatrix[11,9], CPMatrix[11,10], CPMatrix[11,11], CPMatrix[11,12], CPMatrix[11,13], CPMatrix[11,14], CPMatrix[11,15], CPMatrix[11,16], CPMatrix[11,17], CPMatrix[11,18], CPMatrix[11,19], CPMatrix[11,20])
        @printf(f, "\$13\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[12,1], sum(CPMatrix[12,2:end]), CPMatrix[12,2], CPMatrix[12,3], CPMatrix[12,4], CPMatrix[12,5], CPMatrix[12,6], CPMatrix[12,7], CPMatrix[12,8], CPMatrix[12,9], CPMatrix[12,10], CPMatrix[12,11], CPMatrix[12,12], CPMatrix[12,13], CPMatrix[12,14], CPMatrix[12,15], CPMatrix[12,16], CPMatrix[12,17], CPMatrix[12,18], CPMatrix[12,19], CPMatrix[12,20])
        @printf(f, "\$14\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[13,1], sum(CPMatrix[13,2:end]), CPMatrix[13,2], CPMatrix[13,3], CPMatrix[13,4], CPMatrix[13,5], CPMatrix[13,6], CPMatrix[13,7], CPMatrix[13,8], CPMatrix[13,9], CPMatrix[13,10], CPMatrix[13,11], CPMatrix[13,12], CPMatrix[13,13], CPMatrix[13,14], CPMatrix[13,15], CPMatrix[13,16], CPMatrix[13,17], CPMatrix[13,18], CPMatrix[13,19], CPMatrix[13,20])
        @printf(f, "\$15\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[14,1], sum(CPMatrix[14,2:end]), CPMatrix[14,2], CPMatrix[14,3], CPMatrix[14,4], CPMatrix[14,5], CPMatrix[14,6], CPMatrix[14,7], CPMatrix[14,8], CPMatrix[14,9], CPMatrix[14,10], CPMatrix[14,11], CPMatrix[14,12], CPMatrix[14,13], CPMatrix[14,14], CPMatrix[14,15], CPMatrix[14,16], CPMatrix[14,17], CPMatrix[14,18], CPMatrix[14,19], CPMatrix[14,20])
        @printf(f, "\$16\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[15,1], sum(CPMatrix[15,2:end]), CPMatrix[15,2], CPMatrix[15,3], CPMatrix[15,4], CPMatrix[15,5], CPMatrix[15,6], CPMatrix[15,7], CPMatrix[15,8], CPMatrix[15,9], CPMatrix[15,10], CPMatrix[15,11], CPMatrix[15,12], CPMatrix[15,13], CPMatrix[15,14], CPMatrix[15,15], CPMatrix[15,16], CPMatrix[15,17], CPMatrix[15,18], CPMatrix[15,19], CPMatrix[15,20])
        @printf(f, "\$17\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[16,1], sum(CPMatrix[16,2:end]), CPMatrix[16,2], CPMatrix[16,3], CPMatrix[16,4], CPMatrix[16,5], CPMatrix[16,6], CPMatrix[16,7], CPMatrix[16,8], CPMatrix[16,9], CPMatrix[16,10], CPMatrix[16,11], CPMatrix[16,12], CPMatrix[16,13], CPMatrix[16,14], CPMatrix[16,15], CPMatrix[16,16], CPMatrix[16,17], CPMatrix[16,18], CPMatrix[16,19], CPMatrix[16,20])
        @printf(f, "\$18\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[17,1], sum(CPMatrix[17,2:end]), CPMatrix[17,2], CPMatrix[17,3], CPMatrix[17,4], CPMatrix[17,5], CPMatrix[17,6], CPMatrix[17,7], CPMatrix[17,8], CPMatrix[17,9], CPMatrix[17,10], CPMatrix[17,11], CPMatrix[17,12], CPMatrix[17,13], CPMatrix[17,14], CPMatrix[17,15], CPMatrix[17,16], CPMatrix[17,17], CPMatrix[17,18], CPMatrix[17,19], CPMatrix[17,20])
        @printf(f, "\$19\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[18,1], sum(CPMatrix[18,2:end]), CPMatrix[18,2], CPMatrix[18,3], CPMatrix[18,4], CPMatrix[18,5], CPMatrix[18,6], CPMatrix[18,7], CPMatrix[18,8], CPMatrix[18,9], CPMatrix[18,10], CPMatrix[18,11], CPMatrix[18,12], CPMatrix[18,13], CPMatrix[18,14], CPMatrix[18,15], CPMatrix[18,16], CPMatrix[18,17], CPMatrix[18,18], CPMatrix[18,19], CPMatrix[18,20])
        @printf(f, "\$20\$ & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f & %.2f   \\\\ \n", CPMatrix[19,1], sum(CPMatrix[19,2:end]), CPMatrix[19,2], CPMatrix[19,3], CPMatrix[19,4], CPMatrix[19,5], CPMatrix[19,6], CPMatrix[19,7], CPMatrix[19,8], CPMatrix[19,9], CPMatrix[19,10], CPMatrix[19,11], CPMatrix[19,12], CPMatrix[19,13], CPMatrix[19,14], CPMatrix[19,15], CPMatrix[19,16], CPMatrix[19,17], CPMatrix[19,18], CPMatrix[19,19], CPMatrix[19,20])
        @printf(f, "\\hline\\hline \n")
        @printf(f, "\\end{tabular} \n \\egroup \n")
        @printf(f, "\\end{table} \n")
        @printf(f, "\\end{landscape} \n")
    end

nothing

end



