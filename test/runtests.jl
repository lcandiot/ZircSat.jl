using ZircSat
using Test

data_path = "/Users/lcandiot/Developer/ZircSat/data/ZircSat_test.csv"
db = "ig"
sys_in = "wt"
delim = ','
header = true
T_0 = 1050.0

@testset "ZircSat.jl" begin
    @test calculate_liquidus_temperature(; T_0, data_path=data_path, db=db, sys_in=sys_in, delim=delim, header=header) ≈ [1059.8, 1059.8] rtol = 1e-1
end
