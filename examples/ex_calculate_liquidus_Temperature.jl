using ZircSat

data_path = "/Users/lcandiot/Developer/MeltMigration_MagmaticSystems/data/experiments/ZircSat_test.csv"
db = "ig"
sys_in = "wt"
delim = ','
header = true
write_csv = true
T_0 = 1000.0

# Run
T = calculate_liquidus_temperature(; T_0, data_path, db, sys_in, delim, header, write_csv)
