using ZircSat
using Test

data_path_Tliq = "../data/ZircSat_test_Tliq_orig.csv"
data_path_Tsat = "../data/ZircSat_test_Tsat_orig.csv"
db = "ig"
sys_in = "wt"
delim = ','
header = true
write_csv = false
T_0 = 1059.0

@testset "ZircSat.jl" begin
    @test calculate_liquidus_temperature(; T_0, data_path=data_path_Tliq, db, sys_in, delim, header, write_csv) ≈ [1059.8] rtol = 1e-1
    @test calculate_ZircSat_temperature(; data_path=data_path_Tsat) ≈ [809.2] rtol = 1e-1
end
