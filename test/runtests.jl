using ZircSat
using Test

data_path = "../data/ZircSat_test_MarxerUlmer2019.csv"
db        = "ig"
sys_in    = "wt"
delim     = ','
header    = true
T_0       = 1092.0

@testset "ZircSat.jl" begin
    @test calculate_liquidus_temperature(; T_0, data_path, db, sys_in, delim, header) ≈ [1093.6] rtol = 1e-1
    @test calculate_ZircSat_temperature(; data_path) ≈ [698.7] rtol = 1e-1
end
