# Functionality of ZircSat
using MAGEMin_C, Printf, DelimitedFiles, DataFramesMeta, CSV

# Liquidus temperature calculation using MAGEMin
function calculate_liquidus_temperature(; T_0::Float64, data_path::String, db::String, sys_in::String, delim::Char, header::Bool, write_csv::Bool)

    # Read in the CSV file
    data, header = readdlm(data_path, delim, header=header)
    df = DataFrame(data, vec(header))

    # Intitialize
    # db         = "ig"                                   # Database selection: ig = igneous (Holland et al., 2018)
    Xoxides    = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"]
    X          = Vector{Float64}(undef, length(Xoxides))
    T_liquidus = Vector{Float64}(undef, length(df[:,1]))

    # Solver opts
    liq_f_max  = 100.0                                  # Melt fraction at liquidus [%]
    max_iter   = 10_000                                 # Maximum no. Newton-Raphson iteration [ ]
    ϵ_tol      = 1e-1                                   # Absolute error tolerance [ ]
    α          = 0.01                                    # Relaxation factor [ ]
    iters      = 0                                      # Iteration counter

    # Loop through data
    for iData in eachindex(df[:,1])
        
        # Load variables, reset, and start MAGEMin
        P    = df[iData, 1    ]                                    # Pressure [kbar]
        X   .= [df[iData, iOx] for iOx in Xoxides] # delta = 0
        T    = copy(T_0)
        data = Initialize_MAGEMin(db, verbose=false, buffer="qfm");

        # Monitor
        println("----------------")
        println("Running data set $iData / $(length(df[:,1]))")
        println(df[iData, :])
        println(" ")

        # Determine liquidus temperature
        err = 1.0; iters = 0
        for iters in 1:max_iter
                    
            # Determine melt fraction at current T
            out = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, B=2.0)

            # Find liq string in MAGEMin phase output and skip iteration if we are still below the solidus
            id_liq = findfirst(x-> x == "liq", out.ph)
            id_fl  = findfirst(x-> x == "fl",  out.ph)
            if isnothing(id_liq)
                # println("Current T < T_solidus. Increasing T by 50 °C")
                T += 10.0
                continue
            elseif out.ph_frac_wt[id_liq] < 0.7
                T += 2.0
            end

            # Calculate misfit
            if isnothing(id_fl)
                liq_f = out.ph_frac_wt[id_liq] * 100.0
            else
                liq_f = (out.ph_frac_wt[id_liq] + out.ph_frac_wt[id_fl]) * 100.0
            end
            
            if liq_f >= 100.0
                T -= 1.0
                println("Liquid Fraction = $(@sprintf("%.8f", liq_f)) %. Correcting T down. New T = $T")
                continue
            end
            err = abs(liq_f_max - liq_f)

            # Stop calculating if error is below tolerance
            if err < ϵ_tol
                break
            end

            # Calculate liquid fraction at perturbed T (necessary for numerical derivative)
            T_p        = T + T * 1e-4
            out_pert   = single_point_minimization(P, T_p, data, X=X, Xoxides=Xoxides, sys_in=sys_in, B=2.0)
            id_liq     = findfirst(x-> x == "liq", out_pert.ph)
            id_fl      = findfirst(x-> x == "fl", out_pert.ph)
            
            if isnothing(id_fl)
                liq_f_pert = out_pert.ph_frac_wt[id_liq] * 100.0
            else
                liq_f_pert = (out_pert.ph_frac_wt[id_liq] + out_pert.ph_frac_wt[id_fl]) * 100.0
            end

            err_pert   = abs(liq_f_max - liq_f_pert)

            # Update T according to Newton-Raphson method
            dErrdT = (err - err_pert) / (T - T_p)
            T     -= α * err / dErrdT

            if iters % 100 == 0
                println("Iteration $(iters):    err = $err; T = $T; liq_f = $liq_f")
            end
        end

        # Print result
        println(" ")
        println("T_liq = $(T)")
        
        # Store liquidus
        T_liquidus[iData] = T

        # MAGEMin final
        Finalize_MAGEMin(data)

    end

    # Add liquidus T to df and save as new .csv
    if write_csv
        insertcols!(df, "T_liq [C]" => T_liquidus)
        data_path_new = replace(data_path, ".csv" => "_Tliq.csv")
        CSV.write(data_path_new, df)
    end

    # Return
    return T_liquidus

end

# Export
export calculate_liquidus_temperature