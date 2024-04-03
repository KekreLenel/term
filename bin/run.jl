include("../src/main.jl")

if ARGS != [] && ARGS[2] == "della"
    const results_folder = "/scratch/gpfs/mlenel/results/"
elseif ARGS != [] && ARGS[2] == "adroit"
    const results_folder = "/scratch/network/mlenel/results/"
else 
    const results_folder = "../output/results/"
end

if ARGS == [] || ~(ARGS[1] == "output")
    
    # rm/create output folders 
    try; rm(results_folder, recursive=true); catch; end;     
    try; mkpath(results_folder);  catch; end;

    main()
end

if ARGS == [] || ARGS[1] == "output"
    makeOutput()
end 

show(to)
