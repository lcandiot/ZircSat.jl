using ZircSat

# Define data
data_path = "./data/ZircSat_test_Tsat_orig.csv"

@show T_sat = calculate_ZircSat_temperature(; data_path)
