# Functionality of ZircSat
using MAGEMin_C, Printf, DelimitedFiles, DataFramesMeta, CSV

# Liquidus temperature calculation using MAGEMin
function calculate_liquidus_temperature(; T_0::Float64, data_path::String, db::String, sys_in::String, delim::Char, header::Bool)

    # Read in the CSV file
    data, header = readdlm(data_path, delim, header=header)
    df = DataFrame(data, vec(header))

    # Check if FeO is total iron
    idFe2O3 = findfirst(x -> x == "Fe2O3", names(df))

    # Intitialize
    if isnothing(idFe2O3)   # Check if FeO is total iron
        data = Initialize_MAGEMin(db, verbose=false, buffer="qfm")
        Xoxides    = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"]
    else
        data = Initialize_MAGEMin(db, verbose=false)
        Xoxides    = ["SiO2"; "Al2O3"; "CaO"; "MgO"; "FeO"; "Fe2O3"; "K2O"; "Na2O"; "TiO2"; "Cr2O3"; "H2O"]
    end
    X          = Vector{Float64}(undef, length(Xoxides))
    T_liquidus = Vector{Float64}(undef, length(df[:,1]))

    # Solver opts
    liq_f_max  = 100.0                                  # Melt fraction at liquidus [%]
    max_iter   = 100                                 # Maximum no. Newton-Raphson iteration [ ]
    ϵ_tol      = 1e-2                                   # Absolute error tolerance [ ]
    α          = 0.7                                    # Relaxation factor [ ]
    iters      = 0                                      # Iteration counter

    # Loop through data
    for iData in eachindex(df[:,1])
        
        # Load variables, reset, and start MAGEMin
        P     = df[iData, 1    ]                                    # Pressure [kbar]
        X    .= [df[iData, iOx] for iOx in Xoxides] # delta = 0
        T     = copy(T_0)
        T_old = 100.0

        # Monitor
        println("----------------")
        println("Running data set $iData / $(length(df[:,1]))")
        println(df[iData, :])
        println(" ")

        # Determine liquidus temperature
        err = 1.0; iters = 0
        for iters in 1:max_iter

            # Determine melt fraction at current T
            if isnothing(idFe2O3)   # Check if FeO is total iron
                out = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in, B=0.0)
            else
                out = single_point_minimization(P, T, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
            end

            # Find liq string in MAGEMin phase output and skip iteration if we are still below the solidus
            id_liq = findfirst(x-> x == "liq", out.ph)
            id_fl  = findfirst(x-> x == "fl",  out.ph)
            if isnothing(id_liq)
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
            if isnothing(idFe2O3)   # Check if FeO is total iron
                out_pert   = single_point_minimization(P, T_p, data, X=X, Xoxides=Xoxides, sys_in=sys_in, B=2.0)
            else
                out_pert   = single_point_minimization(P, T_p, data, X=X, Xoxides=Xoxides, sys_in=sys_in)
            end
            id_liq     = findfirst(x-> x == "liq", out_pert.ph)
            id_fl      = findfirst(x-> x == "fl", out_pert.ph)
            
            if isnothing(id_fl)
                liq_f_pert = out_pert.ph_frac_wt[id_liq] * 100.0
            else
                liq_f_pert = (out_pert.ph_frac_wt[id_liq] + out_pert.ph_frac_wt[id_fl]) * 100.0
            end

            err_pert   = abs(liq_f_max - liq_f_pert)

            # Update T according to Newton-Raphson method
            T_old  = T
            dErrdT = (err - err_pert) / (T - T_p)
            T     -= α * err / dErrdT

            if iters % 10 == 0
                println("Iteration $(iters):    err = $err; T = $T; T_old = $T_old liq_f = $liq_f")
            end

            # Stop calculating if the temperature value does not significantly change anymore
            if abs(T - T_old) < ϵ_tol
                break
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

    # Return
    return T_liquidus

end

# Zircon saturation temperature
function calculate_ZircSat_temperature(; data_path::String)
    # Define elemental optical basicities
    ElementalOpt = Dict("K2O" => 1.4,
    "Na2O" => 1.15,
    "CaO" => 1.0,
    "MgO" => 0.78,
    "MnO" => 1.0,
    "FeO" => 1.0,
    "Al2O3" => 0.60,
    "SiO2" => 0.48,
    "P2O5" => 0.33)

    MolecularMass = Dict("SiO2"  => 60.084,
                         "K2O"   => 94.196,
                         "CaO"   => 56.077,
                         "MgO"   => 40.304,
                         "MnO"   => 70.937,
                         "FeO"   => 71.844,
                         "Al2O3" => 101.961,
                         "H2O"   => 18.015
    )

    df = CSV.read(data_path, DataFrame)

    # Some sanity checks
    minZr = minimum(df.Zr)
    maxZr = maximum(df.Zr)
    if minZr < 10.0 || maxZr > 1000.0
        error("Zirconium content is out of bounds. Please, double check your data input.")
    end

    # Pressure in GPa
    P    = df[:, "P [kbar]"] ./ 10.0

    # Calculate optical basicity and xH2O
    opt  = calculate_optical_basicity(; df, MolecularMass, ElementalOpt)
    xH2O = calculate_water_molfrac(; df, MolecularMass)
    
    T_sat   = @. -5790.0 / (log10(df.Zr) - 0.96 + 1.28 * P - 12.39 * opt - 0.83 * xH2O - 2.06 * P * opt)
    T_sat .-= 273.15

    println("T_sat = $T_sat")

    # Return
    return T_sat
end

# Optical basicity
function calculate_optical_basicity(; df::DataFrame, MolecularMass::Dict, ElementalOpt::Dict)
    # Convert Fe2O3 to FeO if present
    if hasproperty(df, "Fe2O3")
        if hasproperty(df, "FeO")
            df[!, "FeO"] .+= df[!, "Fe2O3"] ./ 1.1111
        else
            df[!, "FeO"]  .= df[!, "Fe2O3"] ./ 1.1111
        end
    end

    # Calculate molar abundances
    SiO2_mol = df.SiO2 ./ MolecularMass["SiO2"]
    K2O_mol = df.K2O ./ MolecularMass["K2O"]
    CaO_mol = df.CaO ./ MolecularMass["CaO"]
    MgO_mol = df.MgO ./ MolecularMass["MgO"]
    MnO_mol = df.MnO ./ MolecularMass["MnO"]
    FeO_mol = df.FeO ./ MolecularMass["FeO"]
    Al2O3_mol = df.Al2O3 ./ MolecularMass["Al2O3"]
    
    total_mol = (SiO2_mol + K2O_mol + CaO_mol + MgO_mol +
                 MnO_mol + FeO_mol + Al2O3_mol)
    
    # Calculate optical basicity
    term1 = (
        ((SiO2_mol  ./ total_mol) .* 2.0 .* ElementalOpt["SiO2"]) .+
        ((K2O_mol   ./ total_mol) .* 1.0 .* ElementalOpt["K2O"])  .+
        ((CaO_mol   ./ total_mol) .* 1.0 .* ElementalOpt["CaO"])  .+
        ((MgO_mol   ./ total_mol) .* 1.0 .* ElementalOpt["MgO"])  .+
        ((MnO_mol   ./ total_mol) .* 1.0 .* ElementalOpt["MnO"])  .+
        ((FeO_mol   ./ total_mol) .* 1.0 .* ElementalOpt["FeO"])  .+
        ((Al2O3_mol ./ total_mol) .* 3.0 .* ElementalOpt["Al2O3"])
    )
    
    term2 = (
        ((SiO2_mol  ./ total_mol) .* 2.0) .+
        ((K2O_mol   ./ total_mol) .* 1.0) .+
        ((CaO_mol   ./ total_mol) .* 1.0) .+
        ((MgO_mol   ./ total_mol) .* 1.0) .+
        ((MnO_mol   ./ total_mol) .* 1.0) .+
        ((FeO_mol   ./ total_mol) .* 1.0) .+
        ((Al2O3_mol ./ total_mol) .* 3.0)
    )
    
    opt_bas = term1 ./ term2
    
    return opt_bas
end

# Water content based on water activity
function calculate_water_molfrac(; df::DataFrame, MolecularMass::Dict)
    
    # Convert Fe2O3 to FeO if present
    if hasproperty(df, "Fe2O3")
        if hasproperty(df, "FeO")
            df[!, "FeO"] .+= df[!, "Fe2O3"] ./ 1.1111
        else
            df[!, "FeO"]  .= df[!, "Fe2O3"] ./ 1.1111
        end
    end

    # Calculate molar weights
    SiO2_mol  = df.SiO2  ./ MolecularMass["SiO2"]
    K2O_mol   = df.K2O   ./ MolecularMass["K2O"]
    CaO_mol   = df.CaO   ./ MolecularMass["CaO"]
    MgO_mol   = df.MgO   ./ MolecularMass["MgO"]
    MnO_mol   = df.MnO   ./ MolecularMass["MnO"]
    FeO_mol   = df.FeO   ./ MolecularMass["FeO"]
    Al2O3_mol = df.Al2O3 ./ MolecularMass["Al2O3"]
    H2O_mol   = df.H2O   ./ MolecularMass["H2O"]

    total_mol = (SiO2_mol + K2O_mol + CaO_mol + MgO_mol +
                 MnO_mol + FeO_mol + Al2O3_mol + H2O_mol)
    
    xH2O = H2O_mol ./ total_mol
    
    return xH2O
end

# Export
export calculate_liquidus_temperature
export calculate_optical_basicity
export calculate_ZircSat_temperature
export calculate_water_molfrac