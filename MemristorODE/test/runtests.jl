using MemristorODE
using Test
using LinearAlgebra

@testset "MemristorODE.jl" begin
     
    @testset "Window Functions" begin
        # Joglekar
        @test window_joglekar(0.5, 0.0, 1.0, 1) == 1.0 - (0.0)^2
        @test window_joglekar(0.0, 0.0, 1.0, 1) == 0.0
        
        # Biolek
        @test window_biolek(0.5, 0.0, 1.0, 1.0, 1) < 1.0
        @test window_biolek(0.0, 0.0, 1.0, 1.0, 1) == 1.0
        @test window_biolek(1.0, 0.0, 1.0, -1.0, 1) == 1.0
    end
 
    @testset "Threshold Model" begin
        p = ThresholdParams()
        @test resistance(0.0, p) == p.R_off
        @test resistance(1.0, p) == p.R_on
        @test current(1.0, 0.5, p) == 1.0 / resistance(0.5, p)
    end
 
    @testset "VTEAM Model" begin
        p = VTEAMParams()
        @test resistance(p.w_on, p) == p.R_on
        @test resistance(p.w_off, p) == p.R_off
        @test conductance(p.w_on, p) == 1.0 / p.R_on
    end
 
    @testset "Signals" begin
        @test sinusoidal_wave(0.0) == 0.0
        @test triangular_wave(0.0) == 0.0
        @test pulse_train(0.0) == -1.0 # Default V_low
    end
 
    @testset "Crossbar Array" begin
        rows, cols = 4, 4
        xbar = CrossbarArray(rows, cols, R_wire=0.0)
         
        V_in = ones(rows)
        I_ideal = ideal_mvm(xbar.conductances, V_in)
        I_sim = simulate_crossbar_mvm(xbar, V_in, apply_noise=false)
         
        @test I_ideal ≈ I_sim
         
        # Test MNA matrix building
        A = build_crossbar_mna_matrix(rows, cols, xbar.conductances, 1.0)
        @test size(A) == (2*rows*cols, 2*rows*cols)
    end
 
end