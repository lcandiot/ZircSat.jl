using ZircSat
using Test

data_path = "../data/ZircSat_test.csv"
db = "ig"
sys_in = "wt"
delim = ','
header = true
write_csv = false
T_0 = 1055.0

@testset "ZircSat.jl" begin
    @test calculate_liquidus_temperature(; T_0, data_path, db, sys_in, delim, header, write_csv) ≈ [1059.8, 1059.8] rtol = 1e-1
end
