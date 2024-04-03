function makeFiguresPaper(data_empirical::RunEmpirics, results::Vector{Run})

    calibs = getDictionary()

    @unpack y_idxs, σ_β, T_irr = results[calibs["benchmark"]].prm
    dur_highD      = Int(round(results[calibs["low_duration"]].sim.mom[88]))
    
    folder = "../output/figures/"

    # Real forward response

    fig_name    = folder*"df1_real.pdf"
    fig_legends = ["", "", "", "", "", ""]
    fig_params  = FigParams{Float64}()

    makeFiguresFwdResponse(data_empirical.df1_real.series[1:19], data_empirical.df1_real.series_lci[1:19], data_empirical.df1_real.series_uci[1:19], [], [], [], [], [], [], fig_name, "f", fig_legends, [], fig_params)

    # Nominal forward response

    fig_name    = folder*"df1.pdf"
    fig_legends = ["", "", "", "", "", ""]
            
    makeFiguresFwdResponse(data_empirical.df1.series, data_empirical.df1.series_lci, data_empirical.df1.series_uci, [], [], [], [], [], [], fig_name, "f", fig_legends, -2.2, fig_params)

    # Nominal yield response

    fig_name    = folder*"dy1.pdf"
    fig_legends = ["", "", "", "", "", ""]
            
    makeFiguresFwdResponse(data_empirical.dy1.series, data_empirical.dy1.series_lci, data_empirical.dy1.series_uci, [], [], [], [], [], [], fig_name, "y", fig_legends, -1.0, fig_params)

    # Scatter plot of change in real forward rates against change in 1-year yield

    fig_name = folder*"df1_real_20.pdf"

    makeScatterPlots([data_empirical.dyscatter], [data_empirical.dfscatter20], [data_empirical.ymscatter], fig_name, fig_params)

    # Scatter plot of change in dealers equity against change in 1-year yield

    fig_name = folder*"dealer.pdf"

    makeScatterPlots([data_empirical.dyscatter_hf, data_empirical.dyscatter], 
    [data_empirical.avg_ret_hf, data_empirical.avg_ret], 
    [data_empirical.ymscatter_hf, data_empirical.ymscatter], fig_name, fig_params)

    # Response to monetary shock

    idx1        = calibs["benchmark"]
    idx2        = calibs["exogenous_wealth"]

    τ_plot      = 19
    t_max       = 62
    r           = normalize1(results[idx1].irM.r, t_max)
    y1          = normalize1(results[idx1].irM.y[y_idxs[1],:], t_max)
    W           = normalize1(log.(results[idx1].irM.W), t_max)
    fwd         = normalize1(results[idx1].irM.fwd[τ_plot,:], t_max)
    Eret        = normalize1(results[idx1].irM.aE[τ_plot+1,1,:], t_max)
    r_noW       = normalize1(results[idx2].irM.r, t_max)
    y1_noW      = normalize1(results[idx2].irM.y[y_idxs[1],:], t_max)
    W_noW       = normalize1(log.(results[idx2].irM.W), t_max)
    fwd_noW     = normalize1(results[idx2].irM.fwd[τ_plot,:], t_max)
    Eret_noW    = normalize1(results[idx2].irM.aE[τ_plot+1,1,:], t_max)
    plot_ir     = [r, y1, W, fwd, fwd-y1, Eret-y1]
    plot_ir_noW = [r_noW, y1_noW, W_noW, fwd_noW, fwd_noW-y1_noW, Eret_noW-y1_noW]
    plot_x      = 0:60
    factor      = 100/abs(y1[1])

    fig_name    = folder*"ir_M.pdf"

    makeFiguresIR(plot_x, plot_ir*factor, plot_ir_noW*factor, fig_name, 1, 0, 0, 0, fig_params)

    # Response to demand shock

    idx1        = calibs["benchmark"]
    idx2        = calibs["exogenous_wealth"]

    β           = normalize1(results[idx1].irβ.β, t_max)
    y1          = normalize1(results[idx1].irβ.y[y_idxs[1],:], t_max)
    W           = normalize1(log.(results[idx1].irβ.W), t_max)
    fwd         = normalize1(results[idx1].irβ.fwd[τ_plot,:], t_max)
    Eret        = normalize1(results[idx1].irβ.aE[τ_plot+1,1,:], t_max)
    β_noW       = normalize1(results[idx2].irβ.β, t_max)
    y1_noW      = normalize1(results[idx2].irβ.y[y_idxs[1],:], t_max)
    W_noW       = normalize1(log.(results[idx2].irβ.W), t_max)
    fwd_noW     = normalize1(results[idx2].irβ.fwd[τ_plot,:], t_max)
    Eret_noW    = normalize1(results[idx2].irβ.aE[τ_plot+1,1,:], t_max)
    plot_ir     = [β, y1, W, fwd, fwd-y1, Eret-y1]
    plot_ir_noW = [β_noW, y1_noW, W_noW, fwd_noW, fwd_noW-y1_noW, Eret_noW-y1_noW]
    plot_x      = 0:60
    factor      = (σ_β/sqrt(12))*10000/abs(β[1])

    fig_name    = folder*"ir_b.pdf"
    
    makeFiguresIR(plot_x, plot_ir*factor, plot_ir_noW*factor, fig_name, 0, 1, 0, 0, fig_params)

    # Response to short-term rate shock

    idx1        = calibs["benchmark"]
    idx2        = calibs["exogenous_wealth"]

    r           = normalize1(results[idx1].irr.r, t_max)
    y1          = normalize1(results[idx1].irr.y[y_idxs[1],:], t_max)
    W           = normalize1(log.(results[idx1].irr.W), t_max)
    fwd         = normalize1(results[idx1].irr.fwd[τ_plot,:], t_max)
    Eret        = normalize1(results[idx1].irr.aE[τ_plot+1,1,:], t_max)
    r_noW       = normalize1(results[idx2].irr.r, t_max)
    y1_noW      = normalize1(results[idx2].irr.y[y_idxs[1],:], t_max)
    W_noW       = normalize1(log.(results[idx2].irr.W), t_max)
    fwd_noW     = normalize1(results[idx2].irr.fwd[τ_plot,:], t_max)
    Eret_noW    = normalize1(results[idx2].irr.aE[τ_plot+1,1,:], t_max)
    plot_ir     = [r, y1, W, fwd, fwd-y1, Eret-y1]
    plot_ir_noW = [r_noW, y1_noW, W_noW, fwd_noW, fwd_noW-y1_noW, Eret_noW-y1_noW]
    plot_x      = 0:60
    factor      = 100/abs(y1[1])

    fig_name    = folder*"ir_r.pdf"
    
    makeFiguresIR(plot_x, plot_ir*factor, plot_ir_noW*factor, fig_name, 0, 0, 0, 0, fig_params)

    # U-shape

    idx1   = calibs["benchmark"]
    idx2   = calibs["exogenous_wealth"]

    df     = normalize2(results[idx1].irM.fwd)
    dfnoW  = normalize2(results[idx2].irM.fwd)
    dy1    = normalize1(results[idx1].irM.y[y_idxs[1],:], 2)
    dy1noW = normalize1(results[idx2].irM.y[y_idxs[1],:], 2)
    dEy1   = results[idx1].irM.aE[1,2:end,2] .- results[idx1].irM.aE[1,2:end,1]

    fig_name    = folder*"df1_real_model_data_M.pdf"
    fig_legends = [L"\text{Data}", L"\text{Model}", L"\xi \rightarrow \infty", L"\Delta E_t y_{t+\tau-1}^{(1)} / \Delta y_t^{(1)}", L"\text{Duration = %$(dur_highD)}", L"\alpha = 0"]

    makeFiguresFwdResponse(data_empirical.df1_real.series[1:19], data_empirical.df1_real.series_lci[1:19], data_empirical.df1_real.series_uci[1:19], df[1:19]./dy1, dfnoW[1:19]./dy1noW, dEy1[1:19]./dy1, [], [], [], fig_name, "f", fig_legends, -0.1, fig_params)

    # U-shape with high persistence of short-term rate

    idx1      = calibs["high_kappa"]

    df_highK  = normalize2(results[idx1].irM.fwd)
    dy1_highK = normalize1(results[idx1].irM.y[y_idxs[1],:], 2)

    fig_name    = folder*"df1_model_M_highK.pdf"
    fig_legends = [L"\text{Data}", L"\text{Model}", L"\xi \rightarrow \infty", L"\Delta E_t y_{t+\tau-1}^{(1)} / \Delta y_t^{(1)}", L"\text{Duration = %$(dur_highD)}", L"\alpha = 0"]

    makeFiguresFwdResponse([], [], [], [], [], [], [], [], df_highK[1:19]./dy1_highK, fig_name, "f", fig_legends, -0.1, fig_params)

    # Response of carry returns to monetary shock

    idx1    = calibs["benchmark"]

    dy1     = normalize1(results[idx1].irM.y[y_idxs[1],:], 2)
    cp5     = results[idx1].irM.CPE[2][4,2:20] - results[idx1].irM.CPE[1][4,2:20]
    cp10    = results[idx1].irM.CPE[2][9,2:20] - results[idx1].irM.CPE[1][9,2:20]
    cp15    = results[idx1].irM.CPE[2][14,2:20] - results[idx1].irM.CPE[1][14,2:20]
    cp20    = results[idx1].irM.CPE[2][19,2:20] - results[idx1].irM.CPE[1][19,2:20]
    plot_ir = [cp5./dy1, cp10./dy1, cp15./dy1, cp20./dy1]
    plot_x  = 1:19

    fig_name = folder*"cp_M_5_20.pdf"

    makeFiguresCP(plot_x, plot_ir, fig_name, fig_params)

    # Figure r star sequence, 5-yr fwd 5-yr term premium

    idx1    = calibs["benchmark"]

    r       = normalize1(results[idx1].irrstar.r, T_irr)
    W       = normalize1(log.(results[idx1].irrstar.W), T_irr)
    tp5f5   = normalize1(results[idx1].irrstar.tp5f5, T_irr)
    plot_ir = [r, W, tp5f5]
    plot_x  = 0:155
    factor  = 10000

    fig_name = folder*"ir_rstar_sequence_tp5f5_w.pdf"

    makeFiguresIR(plot_x, plot_ir*factor, [], fig_name, 0, 0, 0, 0, fig_params)

    # U-shape under alternative calibrations

    idx1     = calibs["benchmark"]
    idx2     = calibs["low_duration"]
    idx3     = calibs["zero_alpha"]

    df       = normalize2(results[idx1].irM.fwd)
    dfhighD  = normalize2(results[idx2].irM.fwd)
    dflowA   = normalize2(results[idx3].irM.fwd)
    dy1      = normalize1(results[idx1].irM.y[y_idxs[1],:], 2)
    dy1highD = normalize1(results[idx2].irM.y[y_idxs[1],:], 2)
    dy1lowA  = normalize1(results[idx3].irM.y[y_idxs[1],:], 2)

    fig_name    = folder*"df1_model_M_alternatives.pdf"
    fig_legends = [L"\text{Data}", L"\text{Model}", L"\xi \rightarrow \infty", L"\Delta E_t y_{t+\tau-1}^{(1)} / \Delta y_t^{(1)}", L"\text{Duration = %$(dur_highD)}", L"\alpha = 0"]

    makeFiguresFwdResponse([], [], [], df[1:19]./dy1, [], [], dfhighD[1:19]./dy1highD, dflowA[1:19]./dy1lowA, [], fig_name, "f", fig_legends, -0.1, fig_params)

    # Response to QE shock

    idx1         = calibs["benchmark"]
    
    β            = normalize1(log.(results[idx1].irQE.β), t_max)
    y10          = normalize1(results[idx1].irQE.y[y_idxs[10],:], t_max)
    y20          = normalize1(results[idx1].irQE.y[y_idxs[20],:], t_max)
    W            = normalize1(log.(results[idx1].irQE.W), t_max)
    fwd10        = normalize1(results[idx1].irQE.fwd[9,:], t_max)
    fwd20        = normalize1(results[idx1].irQE.fwd[19,:], t_max)
    β_lowW       = normalize1(log.(results[idx1].irQE_lowW.β), t_max)
    y10_lowW     = normalize1(results[idx1].irQE_lowW.y[y_idxs[10],:], t_max)
    y20_lowW     = normalize1(results[idx1].irQE_lowW.y[y_idxs[20],:], t_max)
    W_lowW       = normalize1(log.(results[idx1].irQE_lowW.W), t_max)
    fwd10_lowW   = normalize1(results[idx1].irQE_lowW.fwd[9,:], t_max)
    fwd20_lowW   = normalize1(results[idx1].irQE_lowW.fwd[19,:], t_max)
    plot_ir      = [β, y10, y20, W, fwd10, fwd20]
    plot_ir_lowW = [β_lowW, y10_lowW, y20_lowW, W_lowW, fwd10_lowW, fwd20_lowW]
    plot_x       = 0:60
    factor       = 10000

    fig_name     = folder*"ir_qe.pdf"
    
    makeFiguresIR(plot_x, plot_ir_lowW*factor, plot_ir*factor, fig_name, 0, 0, 1, 0, fig_params)

    # Short-term rate vs demand shocks

    idx1 = calibs["benchmark"]

    dfr  = normalize2(results[idx1].irM.fwd)
    dfβ  = normalize2(results[idx1].irβ.fwd)
    dy1r = normalize1(results[idx1].irM.y[y_idxs[1],:], 2)
    dy1β = normalize1(results[idx1].irβ.y[y_idxs[1],:], 2)

    fig_name = folder*"df1_model_M_demand.pdf"
    
    makeFiguresNonFOMC(dfr./dy1r, dfβ./dy1β, fig_name, fig_params)

    # Figure TIPS issues

    fig_name = folder*"tips_issues.pdf"

    makeFiguresTIPS(data_empirical.tips_issues, fig_name, fig_params)

    # Figure Duration and Term Premia

    ploty1 = replace(data_empirical.tp_5yf5y, missing => NaN)
    ploty2 = replace(data_empirical.log_dur_pd, missing => NaN)
    ploty3 = replace(data_empirical.aggincgap, missing => NaN)

    ploty  = [ploty1, ploty2, ploty3]
    xlabs  = data_empirical.yduration

    fig_name = folder*"duration_and_term.pdf"

    makeFiguresDuration(xlabs, ploty, fig_name, fig_params)

    # Figure 10Y Equivalent

    fig_name = folder*"10ye.pdf"

    makeFigures10YE(data_empirical.qe1_10ye, fig_name, fig_params)

    nothing
end

function makeTablesPaper(results::Vector{Run})

    nprms  = length(results)
    calibs = getDictionary()
    
    folder = "../output/tables/"

    idx1 = calibs["benchmark"]
    # calibration
    filename = folder*"calibration.tex"
    makeCalibTable(results[idx1].sim.mom, results[idx1].sim.FB, results[idx1].irQE_lowW, results[idx1].prm, filename)
    
    # state dependence
    filename = folder*"duration_dep.tex"
    makeDurTable(results[idx1].irM_sdreg, filename)

    # long yield moments
    idx2 = calibs["benchmark_exogenous_wealth"]
    filename = folder*"vols.tex"
    makeLongYieldsTable(results[idx1].sim.mom,results[idx2].sim.mom, filename)
    
    # comparison to analytical moments
    idx3 = calibs["exogenous_wealth"]
    filename = folder*"analytical.tex"
    makeAnalyticalTable(results[idx3].sim.mom, results[idx3].sim.FB, results[idx3].prm, filename)
    
    nothing

end

function makeAddlTables(results::Vector{Run})

    nprms  = length(results)
    calibs = getDictionary()
    
    folder = "../output/tables/addl_tables/"
    for ppp = 1:nprms

        calib = string(ppp)

        filename = folder*"moments_"*calib*".tex"
        makeCalibTable(results[ppp].sim.mom, results[ppp].sim.FB, results[ppp].irQE_lowW, results[ppp].prm, filename)

        # numerical parameters
        filename = folder*"numerical_settings_"*calib*".tex"
        makeNumericalTable(results[ppp].prm, filename)

        # additional moments 
        filename = folder*"addl_"*calib*".tex"
        makeAddlTable(results[ppp], results[ppp].prm, filename)

        # regression moments
        filename = folder*"yieldmoms_"*calib*".tex"
        makeYieldMomsTable(results[ppp].sim.mom, filename)
#        
        # analyical solution for exogenous wealth runs
        if results[ppp].prm.endo_wealth == 0
            filename = folder*"analytical_moments_"*calib*".tex"
            makeAnalyticalTable(results[ppp].sim.mom, results[ppp].sim.FB, results[ppp].prm, filename)
        end

    end

    nothing

end



function makeAddlFigures(results::Vector{Run})

    nprms = length(results)
    fig_params  = FigParams{Float64}()
    
    folder = "../output/figures/addl_figures/"

    for nnn = 1:nprms

        calib = string(nnn) 
        @unpack endo_wealth, τ_max, τ, Δτ, y_idxs, σ_β, T_irr, Wbnds, rbnds, βbnds = results[nnn].prm
    
        # Response to monetary shock

        τ_plot  = 19
        t_max   = 62
        state   = normalize1(results[nnn].irM.r, t_max)
        y1      = normalize1(results[nnn].irM.y[y_idxs[1],:], t_max)
        W       = normalize1(log.(results[nnn].irM.W), t_max)
        fwd     = normalize1(results[nnn].irM.fwd[τ_plot,:], t_max)
        Eret    = normalize1(results[nnn].irM.aE[τ_plot+1,1,:], t_max)
        plot_ir = [state, y1, W, fwd, fwd-y1, Eret-y1]
        plot_x  = 0:60
        factor  = 100/abs(y1[1])

        fig_name   = folder*"IR_M_"*calib*".pdf"
        fig_titles = [L"r", L"y^{(1)}", L"Log(W)", L"f^{(19,20)}", L"f^{(19,20)}-y^{(1)}", L"E[r^{(20)}]-y^{(1)}"]

        makeAddlFiguresIR(plot_x, plot_ir*factor, fig_name, fig_titles, 0, fig_params)

        # Response to short-term rate shock

        state   = normalize1(results[nnn].irr.r, t_max)
        y1      = normalize1(results[nnn].irr.y[y_idxs[1],:], t_max)
        W       = normalize1(log.(results[nnn].irr.W), t_max)
        fwd     = normalize1(results[nnn].irr.fwd[τ_plot,:], t_max)
        Eret    = normalize1(results[nnn].irr.aE[τ_plot+1,1,:], t_max)
        plot_ir = [state, y1, W, fwd, fwd-y1, Eret-y1]
        plot_x  = 0:60
        factor  = 100/abs(y1[1])

        fig_name   = folder*"IR_r_"*calib*".pdf"
        fig_titles = [L"r", L"y^{(1)}", L"Log(W)", L"f^{(19,20)}", L"f^{(19,20)}-y^{(1)}", L"E[r^{(20)}]-y^{(1)}"]

        makeAddlFiguresIR(plot_x, plot_ir*factor, fig_name, fig_titles, 0, fig_params)

        # Response to demand shock

        state   = normalize1(results[nnn].irβ.β, t_max)
        y1      = normalize1(results[nnn].irβ.y[y_idxs[1],:], t_max)
        W       = normalize1(log.(results[nnn].irβ.W), t_max)
        fwd     = normalize1(results[nnn].irβ.fwd[τ_plot,:], t_max)
        Eret    = normalize1(results[nnn].irβ.aE[τ_plot+1,1,:], t_max)
        plot_ir = [state, y1, W, fwd, fwd-y1, Eret-y1]
        plot_x  = 0:60
        factor  =  10000

        fig_name   = folder*"IR_b_"*calib*".pdf"
        fig_titles = [L"\beta", L"y^{(1)}", L"Log(W)", L"f^{(19,20)}", L"f^{(19,20)}-y^{(1)}", L"E[r^{(20)}]-y^{(1)}"]

        makeAddlFiguresIR(plot_x, plot_ir*factor, fig_name, fig_titles, 0, fig_params)

        # Response to pure wealth shock

        state   = normalize1(log.(results[nnn].irW.W), t_max)
        y1      = normalize1(results[nnn].irW.y[y_idxs[1],:], t_max)
        W       = normalize1(log.(results[nnn].irW.W), t_max)
        fwd     = normalize1(results[nnn].irW.fwd[τ_plot,:], t_max)
        Eret    = normalize1(results[nnn].irW.aE[τ_plot+1,1,:], t_max)
        plot_ir = [state, y1, W, fwd, fwd-y1, Eret-y1]
        plot_x  = 0:60
        factor  = 1000/abs(state[1])

        fig_name   = folder*"IR_W_"*calib*".pdf"
        fig_titles = [L"Log(W)", L"y^{(1)}", L"Log(W)", L"f^{(19,20)}", L"f^{(19,20)}-y^{(1)}", L"E[r^{(20)}]-y^{(1)}"]
        
        makeAddlFiguresIR(plot_x, plot_ir*factor, fig_name, fig_titles, 0, fig_params)

        # Response to r-star shock

        state   = normalize1(results[nnn].irrstar.r, T_irr)
        W       = normalize1(log.(results[nnn].irrstar.W), T_irr)
        fwd     = normalize1(results[nnn].irrstar.fwd[τ_plot,:], T_irr)
        Eret    = normalize1(results[nnn].irrstar.aE[τ_plot+1,1,:], T_irr)
        y1      = normalize1(results[nnn].irrstar.y[y_idxs[1],:], T_irr)
        plot_ir = [state, W, fwd, Eret, fwd-y1, Eret-y1]
        plot_x  = 0:155
        factor  = 10000

        fig_name   = folder*"IR_rstar_"*calib*".pdf"
        fig_titles = [L"r", L"Log(W)", L"f^{(19,20)}", L"E[r^{(20)}]", L"f^{(19,20)}-y^{(1)}", L"E[r^{(20)}]-y^{(1)}"]

        makeAddlFiguresIR(plot_x, plot_ir*factor, fig_name, fig_titles, 0, fig_params)

        # Response to QE shock

        state   = normalize1(results[nnn].irQE.β, t_max)
        y10     = normalize1(results[nnn].irQE.y[y_idxs[10],:], t_max)
        y20     = normalize1(results[nnn].irQE.y[y_idxs[20],:], t_max)
        W       = normalize1(log.(results[nnn].irQE.W), t_max)
        fwd10   = normalize1(results[nnn].irQE.fwd[9,:], t_max)
        fwd20   = normalize1(results[nnn].irQE.fwd[19,:], t_max)
        plot_ir = [state, y10, y20, W, fwd10, fwd20]
        plot_x  = 0:60
        factor  = 10000

        fig_name   = folder*"IR_QE_"*calib*".pdf"
        fig_titles = [L"\beta", L"y^{(10)}", L"y^{(20)}", L"Log(W)", L"f^{(9,10)}", L"f^{(19,20)}"]

        makeAddlFiguresIR(plot_x, plot_ir*factor, fig_name, fig_titles, 0, fig_params)

        # Response to QE shock starting from low wealth
        try 
        state   = normalize1(results[nnn].irQE_lowW.β, t_max)
        y10     = normalize1(results[nnn].irQE_lowW.y[y_idxs[10],:], t_max)
        y20     = normalize1(results[nnn].irQE_lowW.y[y_idxs[20],:], t_max)
        W       = normalize1(log.(results[nnn].irQE_lowW.W), t_max)
        fwd10   = normalize1(results[nnn].irQE_lowW.fwd[9,:], t_max)
        fwd20   = normalize1(results[nnn].irQE_lowW.fwd[19,:], t_max)
        plot_ir = [state, y10, y20, W, fwd10, fwd20]
        plot_x  = 0:60
        factor  = 10000

        fig_name   = folder*"IR_QE_lowW_"*calib*".pdf"
        fig_titles = [L"\beta", L"y^{(10)}", L"y^{(20)}", L"Log(W)", L"f^{(9,10)}", L"f^{(19,20)}"]

        makeAddlFiguresIR(plot_x, plot_ir*factor, fig_name, fig_titles, 0, fig_params)
        catch me 
        end
        
        # Forward responses to each shock

        dfwd_r  = normalize2(results[nnn].irr.fwd)
        dfwd_β  = normalize2(results[nnn].irβ.fwd)
        dfwd_W  = normalize2(results[nnn].irW.fwd)
        dfwd_M  = normalize2(results[nnn].irM.fwd)
        dy1_r   = normalize1(results[nnn].irr.y[y_idxs[1],:], 2)
        dy1_M   = normalize1(results[nnn].irM.y[y_idxs[1],:], 2)
        dβ      = normalize1(results[nnn].irβ.β, 2)
        dW      = normalize1(log.(results[nnn].irW.W), 2)

        fig_names  = folder.*["df_r_"*calib*".pdf", "df_b_"*calib*".pdf", "df_w_"*calib*".pdf", "df_M_"*calib*".pdf"] 
        fig_titles = [L"df^{(\tau-1,\tau)}/dy^{(1)}", L"df^{(\tau-1,\tau)}/d\beta", L"df^{(\tau-1,\tau)}/dLog(W)", L"df^{(\tau-1,\tau)}/dy^{(1)}"]
        plot_diff  = [dfwd_r./dy1_r, dfwd_β./dβ, dfwd_W./dW, dfwd_M./dy1_M]

        nplots = size(plot_diff,1)
        for ppp = 1:nplots
            makeAddlFiguresTermStructure(plot_diff[ppp], fig_names[ppp], fig_titles[ppp], 0, fig_params)
        end

        # Average yield term structure
        fig_name  = folder*"y_avg_"*calib*".pdf"
        fig_title = L"y^{(\tau)}"

        makeAddlFiguresTermStructure(vec(results[nnn].sim.y_avrg), fig_name, fig_title, 1, fig_params)

        
        if endo_wealth == 1
            # Wealth distribution

            fig_name  = folder*"hist_logW_"*calib*".pdf"
            fig_title = L"\log W"
    
            makeAddlFiguresHistograms(log.(results[nnn].sim.W), fig_name, fig_title, fig_params, Wbnds)

            # Wealth distribution under measure Q

            fig_name  = folder*"Qhist_logW_"*calib*".pdf"
            fig_title = L"\log W"
    
            makeAddlFiguresHistograms(log.(results[nnn].simQ.W), fig_name, fig_title, fig_params, Wbnds)
        end 
        #
        # r distribution 

        fig_name  = folder*"hist_r_"*calib*".pdf"
        fig_title = L"r"
    
        makeAddlFiguresHistograms(results[nnn].sim.r, fig_name, fig_title, fig_params, rbnds)

        # r distribution under measure Q

        fig_name  = folder*"Qhist_r_"*calib*".pdf"
        fig_title = L"r"
    
        makeAddlFiguresHistograms(results[nnn].simQ.r, fig_name, fig_title, fig_params, rbnds)
        #
        # Beta distribution under measure Q

        fig_name  = folder*"hist_b_"*calib*".pdf"
        fig_title = L"\beta"
    
        makeAddlFiguresHistograms(results[nnn].sim.β, fig_name, fig_title, fig_params, βbnds)

        # Beta distribution under measure Q

        fig_name  = folder*"Qhist_b_"*calib*".pdf"
        fig_title = L"\beta"
    
        makeAddlFiguresHistograms(results[nnn].simQ.β, fig_name, fig_title, fig_params, βbnds)

    end

end


