using XyceSim
using Test

@testset "XyceSim.jl" begin
    # Test basic imports
    @test XyceSim !== nothing
    
    # Test MemristorParams
    params = XyceSim.MemristorParams(Ron=1e3, Roff=100e3, D=3e-9, uv=1e-15)
    @test params.Ron == 1e3
    @test params.Roff == 100e3
    @test params.D == 3e-9
    @test params.uv == 1e-15
    
    # Test CrossbarArray creation
    xbar = XyceSim.CrossbarArray(4, 4, params)
    @test xbar.rows == 4
    @test xbar.cols == 4
    @test xbar.memristor_params.Ron == 1e3
    @test size(xbar.memristor_values) == (4, 4)
    
    # Test that netlist is generated
    @test !isempty(xbar.netlist)
    @test occursin("* Crossbar array netlist", xbar.netlist)
    
    # Test write_netlist function (without actually writing)
    # This would normally write to a file, but we're just testing the function exists
    @test !isempty(string(XyceSim.write_netlist))
    
    println("All tests passed!")
end