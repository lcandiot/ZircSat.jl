using ZircSat

# Define data
data_path = "./data/ZircSat_test_Tliq_orig_Tliq.csv"
write_csv = true

@show T_sat = calculate_ZircSat_temperature(; data_path, write_csv)
