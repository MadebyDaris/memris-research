# MemristorODE

Pure Julia ODE-based memristor simulation package.

## Features

- Threshold-based memristor model
- VTEAM (industry standard) memristor model  
- Crossbar array simulation with:
  - Ideal matrix-vector multiplication
  - Device variability and read noise
  - IR drop analysis using sparse Modified Nodal Analysis (MNA)
- Window functions (Joglekar, Biolek)
- Voltage waveform generators (triangular, sinusoidal, pulse train)
- Comprehensive test suite

## Installation

```julia
using Pkg
Pkg.add("MemristorODE")
```

## Usage

```julia
using MemristorODE

# Create threshold model parameters
tp = ThresholdParams(
    R_on = 1e3, R_off = 100e3,
    Vth_p = 0.5, Vth_n = -0.5,
    k_on = 25.0, k_off = 25.0,
    w_init = 0.5, p_window = 1
)

# Simulate single memristor
voltage_func(t) = 1.5 * sin(2π * 1.0 * t)  # 1.5V amplitude, 1Hz frequency
res = simulate_memristor(tp, (0.0, 2.0), voltage_func)

# Create and simulate crossbar array
xbar = CrossbarArray(4, 4, R_on=1e3, R_off=100e3)
V_in = [0.2, 0.5, 0.1, 0.3]  # Input voltages
I_out = simulate_crossbar_mvm(xbar, V_in)  # Ideal MVM
I_out_ir = simulate_crossbar_mvm_with_ir(xbar, V_in)  # With IR drop
```

## Documentation

See the `docs/` directory for detailed documentation generation scripts.

## Testing

Run tests with:
```julia
using Pkg
Pkg.test("MemristorODE")
```