# packages 
include("packages.jl")

# no display
ENV["GKSwstype"]="nul"   
 
# function files
include("params.jl")
include("empirical.jl")
include("solution.jl")
include("results.jl")
include("analytical.jl")
include("output.jl")
include("tables.jl")
include("figures.jl")

gr() # TODO: switch everything to CAIRO

# global timer 
const to = TimerOutput()

function main()

    data_empirical = ImportData()
    
    # get parameterizations 
    allprms = getParams()

    if ARGS == []
        runindices = eachindex(allprms)
    else
        runindices = parse(Int64,ARGS[1])
    end

    for ppp = runindices
    
        file_name = results_folder*"results_"*string(ppp)*".jld2"
        try; rm(file_name); catch; end;
        
        println("===================================\n")
        println("Solving calibration "*string(ppp)*".\n")

        @timeit to "solution" begin
        println("Solve for equilibrium prices.")
        pitp, sitp, P, states, drift_vecs, r0, β0, logW0, nquad, weights, quadshocks, states, 
        shocks, drift_r_vec, drift_β_vec, drift_W_vec, Ereturn_vec, η_r_vec, η_β_vec = calcP(allprms[ppp])
        end 

        println("Simulate model under Q.")
        simQ = simSeriesQ(allprms[ppp], sitp, ppp)
        println("Done.\n")

        println("Expected carry returns.")
        E_ret, aE_ret, fwds, carry_ret, aEitp, sigmaitp = calcCarry(allprms[ppp], P, pitp, sitp, states, r0, β0, logW0) 
        println("Done.\n")
        
        println("Find stoch steady state.")
        ss = calcSS(allprms[ppp], pitp, sitp, aEitp, sigmaitp)
        println("Done.\n")
   
        @timeit to "simulation" begin
        println("Simulate model.")
        sim = simSeries(allprms[ppp], pitp, sitp, aEitp, sigmaitp, ss.s[3])
        println("Done.\n")
        end
        println("")

        W_avrg = mean(sim.W)

        @timeit to "IR" begin
        # baseline IRFs
        @unpack r_bar, β_bar = allprms[ppp]
        initvec = [r_bar, β_bar, W_avrg]

        println("Solve for IRF - no shock")
        shocks = zeros(3)
        ir0 = calcIR(allprms[ppp], pitp, sitp, aEitp, sigmaitp, initvec, shocks, sim, 0.0, 1, [])
        println("Done.\n")

        println("Solve for IRF - r")
        shocks    = zeros(3)
        shocks[1] = 0.2 # r shock
        irr       = calcIR(allprms[ppp], pitp, sitp, aEitp, sigmaitp, initvec, shocks, sim, 0.0, 1, ir0)
        println("Done.\n")

        println("Solve for IRF - β")
        shocks    = zeros(3)
        shocks[2] = 1.0 # β shock
        irβ       = calcIR(allprms[ppp], pitp, sitp, aEitp, sigmaitp, initvec, shocks, sim, 0.0, 1, ir0) 
        println("Done.\n")

        println("Solve for IRF - W")    
        shocks    = zeros(3)
        shocks[3] = -0.01*W_avrg # W shock
        irW       = calcIR(allprms[ppp], pitp, sitp, aEitp, sigmaitp, initvec, shocks, sim, 0.0, 1, ir0) 
        println("Done.\n")
        end

        @timeit to "IRM" begin
        @unpack calc_IRM, calc_QE = allprms[ppp]
        if calc_IRM == true || calc_QE == true

            println("Solve for prices with demand and interest rate shifts.")
            println("No shift.")
            pitp0_vec, sitp0_vec, aEitp0_vec, r0_shiftvec, β0_shiftvec = calcPx(allprms[ppp], pitp, sitp, drift_vecs, nquad, weights, quadshocks, [],[], W_avrg)
           
            if calc_IRM == true

                println("r shift.")
                pitpr_vec, sitpr_vec, aEitpr_vec, rr_shiftvec, βr_shiftvec = calcPx(allprms[ppp], pitp, sitp, drift_vecs, nquad, weights, quadshocks, 0.2, [], W_avrg)    
                println("Solve for IRF - rM")
                initvec   = [r_bar, β_bar, W_avrg]
                shocks    = zeros(3)
                ir0       = calcIRM(allprms[ppp], pitp, sitp, aEitp, pitp0_vec, sitp0_vec, aEitp0_vec, r0_shiftvec, β0_shiftvec, initvec, shocks, sim, 0.0, 1, [])
                irM       = calcIRM(allprms[ppp], pitp, sitp, aEitp, pitpr_vec, sitpr_vec, aEitpr_vec, rr_shiftvec, βr_shiftvec, initvec, shocks, sim, 0.0, 1, ir0)
                irM_sdreg = calcRegIRM(allprms[ppp], pitp, pitpr_vec, pitp0_vec, initvec, sim, string(ppp))
                println("Done.\n")

            else

                irM       = ir0
                irM_sdreg = StateDepReg(0, 0, 0, 0, 0)

            end

            if calc_QE == true

                println("QE shift.")
                pitpQE_vec, sitpQE_vec, aEitpQE_vec, rQE_shiftvec, βQE_shiftvec = calcPx(allprms[ppp], pitp, sitp, drift_vecs, nquad, weights, quadshocks, [], data_empirical.qe_purchases, W_avrg)    

                println("Solve for QE")
                initvec = [r_bar, β_bar, W_avrg]
                shocks  = zeros(3)
                ir0     = calcIRM(allprms[ppp], pitp, sitp, aEitp, pitp0_vec, sitp0_vec, aEitp0_vec, r0_shiftvec, β0_shiftvec, initvec, shocks, sim, 0.0, 0, [])
                irQE    = calcIRM(allprms[ppp], pitp, sitp, aEitp, pitpQE_vec, sitpQE_vec, aEitpQE_vec, rQE_shiftvec, βQE_shiftvec, initvec, shocks, sim, 0.0, 0, ir0)

                println("Solve for QE - low W")
                initvec   = [r_bar, β_bar, W_avrg*2.0/3.0]
                shocks    = zeros(3)
                ir0       = calcIRM(allprms[ppp], pitp, sitp, aEitp, pitp0_vec, sitp0_vec, aEitp0_vec, r0_shiftvec, β0_shiftvec, initvec, shocks, sim, 0.0, 0, [])
                irQE_lowW = calcIRM(allprms[ppp], pitp, sitp, aEitp, pitpQE_vec, sitpQE_vec, aEitpQE_vec, rQE_shiftvec, βQE_shiftvec, initvec, shocks, sim, 0.0, 0, ir0)


            else

                irQE      =  ir0
                irQE_lowW =  ir0

            end

        else

            irM       = ir0
            irM_sdreg = StateDepReg(0, 0, 0, 0, 0)
            irQE      =  ir0
            irQE_lowW =  ir0

        end
        end 

        @timeit to "IRR" begin
        @unpack r_bar_shift = allprms[ppp]
        initvec = [r_bar, β_bar, W_avrg]
        ir0     = calcIRR(allprms[ppp], pitp, sitp, aEitp, initvec, 0.0, [], zeros(size(r_bar_shift,1)))
        irrstar = calcIRR(allprms[ppp], pitp, sitp, aEitp, initvec, 0.0, ir0, r_bar_shift)
        println("Done.\n")
        end 
        
        tmp_results = Run(ss, irM_sdreg, sim, simQ, irr, irβ, irW, irM, irQE, irQE_lowW, irrstar, allprms[ppp])
        FileIO.save(file_name,"tmp_results",tmp_results)
        GC.gc()

    end
   
end

function makeOutput() 

    allprms = getParams()
    data_empirical = ImportData()
    
    try; rm("../output/tables", recursive=true); catch; end;     
    try; rm("../output/figures", recursive=true); catch; end;     
    try; mkpath("../output/tables/addl_tables"); catch; end;     
    try; mkpath("../output/figures/addl_figures"); catch; end;     

    default_file = results_folder*"results_1.jld2"

    results = Vector{Run}()
    for ppp in eachindex(allprms)
        
        file_name = results_folder*"results_"*string(ppp)*".jld2"

        if isfile(file_name)
            results_temp = FileIO.load(file_name,"tmp_results")
        else
            println("WARNING: Not all result files present")
            results_temp = FileIO.load(default_file,"tmp_results")
        end

        push!(results, results_temp)
    end
    
    # Additional figures
    makeAddlFigures(results)

    # Additional Tables
    makeAddlTables(results)

    # Figures for paper
    if size(allprms,1) >= 7
        makeTablesPaper(results)
        makeFiguresPaper(data_empirical, results)
    end


end



