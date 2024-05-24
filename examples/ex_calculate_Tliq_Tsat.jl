using ZircSat, CSV, DelimitedFiles, DataFramesMeta

data_path = "/Users/lcandiot/Developer/ZircSat/data/ZircSat_test_MarxerUlmer2019.csv"   # Path to your data
db        = "ig"                                                                        # MAGEMin data base - MUST BE SET TO "ig" FOR CURRENT RELEASE
sys_in    = "wt"                                                                        # Weight percent: "wt", mol percent: "mol". It is recommend to give input as "wt"
delim     = ','                                                                         # Delimiter of your data file
header    = true                                                                        # Does your data file have a header? - It needs one!
write_csv = true                                                                        # Write output into a new CSV file
T_0       = 800.0                                                                      # Initial guess for Newton-Raphson. The closer T_0 to T_liq the faster the convergence

# Run
T_liq = calculate_liquidus_temperature(; T_0, data_path, db, sys_in, delim, header)
T_sat = calculate_ZircSat_temperature(; data_path)

# Save
if write_csv
    # Read in the old CSV data
    data, header = readdlm(data_path, delim, header=header)
    df = DataFrame(data, vec(header))

    # Append columns for temperatures
    insertcols!(df, "T_liq [C]"  => T_liq)
    insertcols!(df, "T_sat [C]"  => T_sat)
    insertcols!(df, "T_diff [C]" => abs.(T_sat .- T_liq))

    # Save new CSV data
    data_path_new = replace(data_path, ".csv" => "_new.csv")
    CSV.write(data_path_new, df)
end
