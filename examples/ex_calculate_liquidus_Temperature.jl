using ZircSat

data_path = "/Users/lcandiot/Developer/ZircSat/data/ZircSat_test_MarxerUlmer2019.csv"
db = "ig"
sys_in = "wt"
delim = ','
header = true
write_csv = true
T_0 = 800.0

# Run
T = calculate_liquidus_temperature(; T_0, data_path, db, sys_in, delim, header)
