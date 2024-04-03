function makeFiguresFwdResponse(data, data_lb, data_ub, model, exog_wealth, expectation_hypothesis, model_duration, model_alpha, model_kappa, fig_name, stringvar, fig_labels, ylb, fig_params)
        
    @unpack_FigParams fig_params 

    if model != []
            ntau = size(model,1)
    elseif model_kappa != []
            ntau = size(model_kappa,1)
    else
            ntau = size(data,1)
    end

    if (model != []) 
            xylim = (1, ntau, ylb, 0.9)
    elseif (model == [] && model_kappa == [] && ylb != [])
            xylim = (1, ntau, ylb, nothing)
    else
            xylim = (1, ntau, nothing, nothing)
    end

    if stringvar == "f"
            fig_title = L"\Delta f^{(\tau-1,\tau)}/ \Delta y^{(1)}"
    else
            fig_title = L"\Delta y^{(\tau)}/ \Delta y^{(1)}"
    end

    fig = Figure()
    ax  = Axis(fig[1,1]; title = fig_title, limits = xylim, xtickalign=1.0, ytickalign=1.0, format1..., axisargs...)
    resize_to_layout!(fig)

    xtick_positions = 4:5:ntau
    xtick_labels = string.(5:5:(ntau+1))
    ax.xticks = (xtick_positions, xtick_labels)

    if data != []

            lines!(ax, 1..ntau, data, label=fig_labels[1], color=:royalblue3, linewidth=linewidth)

            band!(ax, 1:ntau, data_lb, data_ub, color=(:royalblue3, 0.5))

    end

    if model != []
            lines!(ax, 1..ntau, model, label=fig_labels[2], color=:navyblue, linewidth=linewidth)
    end

    if exog_wealth != []
            lines!(ax, 1..ntau, exog_wealth, label=fig_labels[3], color=:navyblue, linestyle=:dot, linewidth=linewidth)
            lines!(ax, 1..ntau, expectation_hypothesis, label=fig_labels[4], color=:black, linewidth=linewidth)
    end

    if model_alpha != []
            lines!(ax, 1..ntau, model_duration, label=fig_labels[5], color=:navyblue, linestyle=:dash, linewidth=linewidth)
            lines!(ax, 1..ntau, model_alpha, label=fig_labels[6], color=:navyblue, linestyle=:dot, linewidth=linewidth)
    end

    if model_kappa != []
            lines!(ax, 1..ntau, model_kappa, color=:navyblue, linewidth=linewidth)
    end

    if model_kappa == []
    hlines!(ax, [0], color=:black, linewidth=0.5)
    end 

    if model_alpha == [] && model_kappa == []
        leg = axislegend(ax; patchsize=(10,4), legargs... )
    elseif model_kappa == []
        leg = axislegend(ax; patchsize=(20,4), legargs... )
    end

    if model == [] && model_kappa == []
            delete!(leg)
    end

    save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeScatterPlots(plotx, ploty, plotannotations, fig_name, fig_params)

    @unpack_FigParams fig_params 

        nsubplots = size(plotx,1)
        fig = Figure()

        if nsubplots == 1
                ylab  = [L"\Delta f^{(19,20)}"]
                ytic  = [-0.2:0.2:0.2]
                xtic  = [-0.2:0.2:0.2]
                xylim = [(-0.3, 0.3, -0.3, 0.3)]
        else
                ylab = ["30-minute change", "One-day change"]
                ytic  = [-4:4:4, -6:11:16]
                xtic  = [-0.2:0.2:0.2, -0.2:0.2:0.2]
                xylim = [(-0.3, 0.3, -4.2, 4.2), (-0.3, 0.3, -6.0, 16.0)]
        end

        for i = 1:Int(nsubplots)

                ax = Axis(fig[Int(ceil(i/2)), Int((i-1)%2 + 1)]; title = "", xtickalign=1.0, ytickalign=1.0,
                xlabel=L"\Delta \hat{y}^{(1)}", ylabel=ylab[i], xticks = xtic[i], yticks = ytic[i], limits = xylim[i], format2..., axisargs... )

                CairoMakie.scatter!(ax, plotx[i], ploty[i], color=:royalblue3, markersize=markersize, strokewidth=strokewidth)

                for iii in 1:size(plotx[i],1)
                        CairoMakie.text!(ax, plotx[i][iii], ploty[i][iii] + 0.01, text=plotannotations[i][iii], fontsize = 8, valign = :bottom)
                end
                hlines!(ax, [0], color=:black, linewidth=0.5)

        end
        resize_to_layout!(fig)
        save(fig_name, fig, pt_per_unit=pt_per_unit)

end


function makeFiguresIR(plotx, ploty1, ploty2, fig_name, M_dummy, b_dummy, QE_dummy, tp10_dummy, fig_params)

    @unpack_FigParams fig_params 
   
    if (size(plotx,1) == 61) & (QE_dummy == 0)

            if b_dummy == 0
                    fig_titles  = [L"r", L"y^{(1)}", L"\log(\!W)", L"f^{(19,20)}", L"f^{(19,20)}-y^{(1)}", L"E[r^{(20)}]-y^{(1)}"]
            else
                    fig_titles  = [L"\beta", L"y^{(1)}", L"\log(\!W)", L"f^{(19,20)}", L"f^{(19,20)}-y^{(1)}", L"E[r^{(20)}]-y^{(1)}"]
            end
            fig_legends = [[L"\text{Model}", L"\xi \rightarrow \infty"], ["", ""], ["", ""], ["", ""], ["", ""], ["", ""]]
            xtick_positions = 10:20:50
            xtick_labels = string.(10:20:50)
            xylim = (0, 60, nothing, nothing)

    elseif (size(plotx,1) == 61) & (QE_dummy == 1)

            fig_titles = [L"\log \int \theta(\tau) d\tau", L"y^{(10)}", L"y^{(20)}", L"\log(\!W)", L"f^{(9,10)}", L"f^{(19,20)}"]
            fig_legends = [[" ", " "], ["", ""], [L"W_0 = 0.67W", L"W_0 = W"], ["", ""], ["", ""], ["", ""]]
            xtick_positions = [0,30,60]
            xtick_labels = ["Apr09", "Oct11", "Apr14   "]
            xylim = (0, 60, nothing, nothing)

    else
            
            if tp10_dummy == 1
                    fig_titles = [L"r", L"\log(\!W)", "10-year TP"]
            else
                    fig_titles = [L"r", L"\log(\!W)", "5-yr fwd, 5-yr TP"]
            end

            fig_legends = [[""], [""], [""]]
            xtick_positions = [0,77,155]
            xtick_labels = ["Jan04", "Jun10", "Dec16  "]
            xylim = (0, 155, nothing, nothing)

    end

    nsubplots = size(fig_titles,1)

    fig = Figure()
    grid = fig[1,1] = GridLayout()

    for i in 1:nsubplots

            ax = Axis(grid[Int(ceil(i / 3)), ((i-1) % 3) + 1]; title = fig_titles[i], limits = xylim, ylabel="bp", xtickalign=1.0, ytickalign=1.0, axisargs..., format4...  )
            lines!(ax, plotx, ploty1[i], label=fig_legends[i][1], color=:navyblue, linewidth=linewidth)

            if ploty2 != []
                    lines!(ax, plotx, ploty2[i], label=fig_legends[i][2], color=:navyblue, linestyle=:dot, linewidth=linewidth)
            end

            hlines!(ax, [0], color=:black, linewidth=0.5)
            ax.xticks = (xtick_positions, xtick_labels)
                            
            if (i == 1) && (size(plotx,1) != 156) && (ploty2 != []) && (QE_dummy == 0)
                    if M_dummy == 1
                            axislegend(ax; patchsize=(8,3), rowgap=2, legargs... )
                    else
                            axislegend(ax; legargs..., position=:rb, patchsize=(8,3), rowgap=2)
                    end
            end

            if (i == 3) && (size(plotx,1) != 156) && (ploty2 != []) && (QE_dummy == 1)
                    axislegend(ax; legargs..., position=:rb, patchsize=(8,3), rowgap=2)
            end

    end

    colgap!(grid,1)
    resize_to_layout!(fig)

    save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeFiguresCP(plot_x, plot_y, fig_name, fig_params)

    @unpack_FigParams fig_params 

        fig_labels = [L"\tau=5", L"\tau=10", L"\tau=15", L"\tau=20"]
        fig_title  = L"\Delta E_{t}\left[r_{t+h}^{(\tau+1-h)}-r_{t+h}^{(\tau-h)}\right]/\Delta y_{t}^{(1)}"
        linstyles  =[:dashdot, :dot, :dash, :solid]

        fig = Figure()

        ax = Axis(fig[1,1]; title = fig_title, limits = (1, 19, nothing, nothing), xtickalign=1.0, ytickalign=1.0, format1..., axisargs...)
        ax.xticks = (4:5:19, string.(5:5:20))

        for i = 1:size(plot_y,1)
                lines!(ax, plot_x, plot_y[i], label=fig_labels[i], color=:navyblue, linewidth=linewidth, linestyle=linstyles[i])
        end

        hlines!(ax, [0], color=:black, linewidth=0.5)
        axislegend(ax; patchsize=(24,6), rowgap=2, legargs...)

    resize_to_layout!(fig)
        save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeFiguresReg(data, data_lb, data_ub, model, no_beta, fig_name, fig_title, fig_params)

    @unpack_FigParams fig_params 
    
        ntau = size(model,1)

        fig = Figure()

        ax = Axis(fig[1,1]; title = fig_title, limits = (1, 19, nothing, nothing), xtickalign=1.0, ytickalign=1.0, format1..., axisargs...)
        ax.xticks = (4:5:ntau, string.(5:5:(ntau+1)))

        lines!(ax, 1..ntau, data, label=L"\text{Data}", color=:royalblue3, linewidth=linewidth)
        band!(ax, 1:ntau, data_lb, data_ub, color=(:royalblue3, 0.5))

        lines!(ax, 1..ntau, model, label=L"\text{Model}", color=:navyblue, linewidth=linewidth)
        lines!(ax, 1..ntau, no_beta, label=L"\sigma_{\beta} = 0", color=:navyblue, linestyle=:dot, linewidth=linewidth)

        hlines!(ax, [0], color=:black, linewidth=0.5)

        axislegend(ax; legargs..., position=:lt, patchsize=(8,8))

    resize_to_layout!(fig)
        save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeFiguresNonFOMC(ploty1, ploty2, fig_name, fig_params)

    @unpack_FigParams fig_params 
        ntau = size(ploty1,1)

        fig = Figure()

        ax = Axis(fig[1,1]; title = L"\Delta f^{(\tau-1,\tau)}/ \Delta y^{(1)}", limits = (1, 19, nothing, nothing), xtickalign=1.0, ytickalign=1.0, axisargs..., format1...)
        ax.xticks = (4:5:ntau, string.(5:5:(ntau+1)))

        lines!(ax, 1..ntau, ploty1, label=L"r", color=:navyblue, linewidth=linewidth)
        lines!(ax, 1..ntau, ploty2, label=L"\beta", color=:navyblue, linestyle=:dash, linewidth=linewidth)

        hlines!(ax, [0], color=:black, linewidth=0.5)

        axislegend(ax; legargs..., position=:lt, patchsize=(24,6))

    resize_to_layout!(fig)
        save(fig_name, fig, pt_per_unit=1.0)

end

function makeFiguresTIPS(ploty, fig_name, fig_params)

    @unpack_FigParams fig_params 
        nt    = size(ploty,1)
        plotx = 1:nt

        fig = Figure()

        ax = Axis(fig[1,1]; title = "Outstanding TIPS", ylabel="Remaining time to maturity", limits = (1, nt, nothing, nothing), xtickalign=1.0, ytickalign=1.0, axisargs..., format1...)
        ax.xticks = (12:36:228, string.(1998:3:2016))
        ax.yticks = (0:5:30, string.(0:5:30))

        for ttt = 1:size(ploty,2)
                ploty_temp = ploty[:,ttt]
                ploty_rest = Float64.(ploty_temp[.~ismissing.(ploty_temp)])
                plotx_rest = Float64.(plotx[.~ismissing.(ploty_temp)])
                lines!(ax, plotx_rest, ploty_rest, color=:royalblue3, linewidth=linewidth)
        end
        hlines!(ax, [0], color=:black, linewidth=0.5)

    resize_to_layout!(fig)
        save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeFigures10YE(ploty, fig_name,fig_params)

    @unpack_FigParams fig_params 
        nt    = size(ploty,1)
        plotx = 1:nt

        fig = Figure()

        ax = Axis(fig[1,1]; title = "Purchases in equivalent 10-year zero coupon bonds",  limits = (1, nt, nothing, nothing),
        ylabel="\$ billion", xtickalign=1.0, ytickalign=1.0, axisargs..., format1...)
        ax.xticks = ([2,5,8,11,14,17], ["Dec08","Mar09","Jun09","Sep09","Dec09","Mar10   "])

        lines!(ax, plotx, ploty, color=:royalblue3, linewidth=linewidth)

        hlines!(ax, [0], color=:black, linewidth=0.5)

    resize_to_layout!(fig)
        save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeFiguresDuration(xlabs, ploty, fig_name,fig_params)

    @unpack_FigParams fig_params 
        nt    = size(xlabs,1)
        plotx = 1:nt

        fig = Figure()

        ax11 = Axis(fig[1,1]; ylabel="%",limits = (1, nt, nothing, nothing), xtickalign=1.0, ytickalign=1.0, axisargs..., format3...)
        ax12 = Axis(fig[1,1]; limits = (1, nt, nothing, nothing), xtickalign=1.0, ytickalign=1.0, yaxisposition = :right, axisargs..., format3...)
        ax21 = Axis(fig[2,1]; limits = (1, nt, nothing, nothing), xtickalign=1.0, ytickalign=1.0, axisargs..., format3...)
        ax22 = Axis(fig[2,1]; limits = (1, nt, nothing, nothing), xtickalign=1, ytickalign=1.0, yaxisposition = :right, axisargs..., format3...)

        ax11.xticks = (1:Int(round((nt-1)/6)):nt, string.(xlabs[1:Int(round((nt-1)/6)):nt]))
        ax12.xticks = (1:Int(round((nt-1)/6)):nt, string.(xlabs[1:Int(round((nt-1)/6)):nt]))
        ax21.xticks = (1:Int(round((nt-1)/6)):nt, string.(xlabs[1:Int(round((nt-1)/6)):nt]))
        ax22.xticks = (1:Int(round((nt-1)/6)):nt, string.(xlabs[1:Int(round((nt-1)/6)):nt]))

        lines!(ax11, plotx[.~isnan.(ploty[2])], ploty[2][.~isnan.(ploty[2])] .- ploty[2][.~isnan.(ploty[2])][1] .+ 1, color=:navyblue, linewidth=linewidth, label="Log dealer dur, left")
        lines!(ax11, plotx, NaN*ploty[1], color=:royalblue3, linewidth=linewidth, label="5-yr fwd, 5-yr TP, right")
        lines!(ax12, plotx[.~isnan.(ploty[1])], ploty[1][.~isnan.(ploty[1])], color=:royalblue3, linewidth=linewidth)

        lines!(ax21, plotx[.~isnan.(ploty[2])], ploty[2][.~isnan.(ploty[2])] .- ploty[2][.~isnan.(ploty[2])][1] .+ 1, color=:navyblue, linewidth=linewidth, label="Log dealer dur, left")
        lines!(ax21, plotx, NaN*ploty[3], color=:royalblue3, linewidth=linewidth, label="Dealer income gap, right")
        lines!(ax22, plotx[.~isnan.(ploty[3])], ploty[3][.~isnan.(ploty[3])], color=:royalblue3, linewidth=linewidth)

        hlines!(ax11, [0], color=:black, linewidth=0.5)
        hlines!(ax21, [0], color=:black, linewidth=0.5)

        axislegend(ax11; legargs..., position=:lt, patchsize = (15.0f0, 10.0f0), rowgap=0.0, padding=(2.0f0, 2.0f0, 2.0f0, 2.0f0))
        axislegend(ax21; legargs..., position=:lt, patchsize = (15.0f0, 10.0f0), rowgap=0.0, padding=(2.0f0, 2.0f0, 2.0f0, 2.0f0))

    resize_to_layout!(fig)
        save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeAddlFiguresIR(plotx, ploty, fig_name, fig_titles, rstar_dummy, fig_params)

    @unpack_FigParams fig_params 

    if rstar_dummy == 0
            xtick_positions = 10:20:50
            xtick_labels = string.(10:20:50)
            xylim = (0, 60, nothing, nothing)
    else
            xtick_positions = [0,77,155]
            xtick_labels = ["Jan04", "Jun10", "Dec16  "]
            xylim = (0, 155, nothing, nothing)
    end

    nsubplots = size(fig_titles,1)

    fig = Figure()
    grid = fig[1,1] = GridLayout()

    for i in 1:nsubplots

            ax = Axis(grid[Int(ceil(i / 3)), ((i-1) % 3) + 1]; title = fig_titles[i], limits = xylim, ylabel="bp", xtickalign=1.0, ytickalign=1.0, axisargs..., format4...)
            lines!(ax, plotx, ploty[i], color=:navyblue, linewidth=linewidth)

            hlines!(ax, [0], color=:black, linewidth=0.5)
            ax.xticks = (xtick_positions, xtick_labels)

    end

    colgap!(grid,1)

    resize_to_layout!(fig)

    save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeAddlFiguresTermStructure(model, fig_name, fig_title, start_from_one, fig_params)

    @unpack_FigParams fig_params 
    
    ntau      = size(model,1)
    xylim     = (1, ntau, nothing, nothing)

    fig = Figure()

    ax = Axis(fig[1,1]; title = fig_title, limits = xylim, xtickalign=1.0, ytickalign=1.0, axisargs..., format1...)

    if start_from_one == 1
            xtick_positions = 5:5:ntau
            xtick_labels = string.(5:5:ntau)
    else
            xtick_positions = 4:5:ntau
            xtick_labels = string.(5:5:(ntau+1))
    end
    ax.xticks = (xtick_positions, xtick_labels)

    lines!(ax, 1..ntau, model, color=:navyblue, linewidth=linewidth)

    hlines!(ax, [0], color=:black, linewidth=0.5)

    resize_to_layout!(fig)

    save(fig_name, fig, pt_per_unit=pt_per_unit)

end

function makeAddlFiguresHistograms(distribution, fig_name, fig_title, fig_params, bounds = [])

    @unpack_FigParams fig_params 

    fig = Figure()

    ax = Axis(fig[1,1]; title = fig_title, xtickalign=1.0, ytickalign=1.0, axisargs..., format1...)

    hist!(ax, distribution, bins = 100, color = :navyblue, strokewidth = 1, strokecolor = :black)
    
    if bounds != []
        vlines!(ax, bounds, color = :red)
    end 

    vlines!(ax, bounds, color = :red)

    hlines!(ax, [0], color=:black, linewidth=0.5)

    resize_to_layout!(fig)

    save(fig_name, fig, pt_per_unit=1.0)

end
