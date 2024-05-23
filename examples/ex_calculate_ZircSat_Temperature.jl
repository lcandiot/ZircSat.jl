using ZircSat

# Define data
data_path = "./data/SampleWR_Data.csv"

@show T_sat = calculate_ZircSat_temperature(; data_path)
