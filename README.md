# Memristor Research: State-of-the-Art Analysis

[![Paper](https://img.shields.io/badge/Paper-PDF-red)](memristor-sota.pdf)

Analysis of memristors and resistive switching for neuromorphic computing and in-memory computation.

## Paper

The full paper (`memristor-sota.pdf`) covers:
- Historical evolution of memristors (Chua's theory to HP Labs 2008)
- Resistive switching mechanisms (SET/RESET, bipolar vs unipolar)
- Crossbar architectures (0T1R passive vs 1T1R active arrays)
- Matrix-Vector Multiplication (MVM) using Ohm's Law + KCL
- IR drop analysis and scalability challenges
- Neuromorphic computing applications

## Repository Structure

```
memris-research/
├── memristor-sota.pdf      # Full paper
├── README.md
├── spice/                  # LTspice simulation files
│   ├── memristor.sub       # Biolek memristor subcircuit model
│   ├── memristor.asy       # LTspice symbol
│   ├── memristorcharac.asc # I-V characterization circuit
│   ├── PassiveCrossbar.asc # 3×3 passive (0T1R) array
│   └── PassiveCrossbar1R1T.asc # 3×3 active (1T1R) array
├── figures/                # Generated analysis plots
│   ├── model_comparison.png
│   ├── ir_drop_map.png
│   ├── scalability_plot.png
│   └── differentiable_sim.png
└── demo.ipynb              # Interactive demonstration notebook
```

## SPICE Models

### Biolek Memristor Model

The memristor subcircuit (`spice/memristor.sub`) implements the Biolek window function model:

```spice
.SUBCKT memristor Plus Minus PARAMS:
+ Ron=1K Roff=100K Rinit=80K D=10N uv=10F p=1

* State variable ODE
Gx 0 x value={ I(Emem)*uv*Ron/D^2*f(V(x),p) }
Cx x 0 1 IC={(Roff-Rinit)/(Roff-Ron)}

* Joglekar window function
.func f(x,p)={1-(2*x-1)^(2*p)}

.ENDS memristor
```

**Parameters:**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `Ron` | Low resistance state (LRS) | 1 kΩ |
| `Roff` | High resistance state (HRS) | 100 kΩ |
| `Rinit` | Initial resistance | 80 kΩ |
| `D` | Thin film width | 10 nm |
| `uv` | Ion mobility coefficient | 10 fm²/Vs |
| `p` | Window function exponent | 1 |

### Running Simulations

1. Open LTspice
2. Add `memristor.sub` to your library path
3. Open `PassiveCrossbar.asc` for passive array simulation
4. Run transient analysis (`.tran 0 2m 0 1u`)

## Key Results

### I-V Hysteresis
The memristor exhibits the characteristic "pinched hysteresis loop" - current passes through origin regardless of voltage polarity.

![Model Comparison](figures/model_comparison.png)

### IR Drop Analysis
Wire resistance causes voltage attenuation in large arrays, limiting practical crossbar size:

![IR Drop Map](figures/ir_drop_map.png)

| Wire Resistance | MVM Error |
|-----------------|-----------|
| 0 Ω | 0% (ideal) |
| 1 Ω | ~2% |
| 10 Ω | ~8% |
| 50 Ω | ~20% |

### Scalability
Simulation timing for crossbar arrays of increasing size:

![Scalability](figures/scalability_plot.png)

### Julia Simulation Framework (WIP)

> **Note:** The Julia-based memristor simulation framework is still under development.

The framework will include:
- `ThresholdParams` / `VTEAMParams` - Device parameter structures
- `simulate_memristor()` - Single device transient simulation
- `CrossbarArray` - Crossbar array with conductance matrix
- `simulate_crossbar_mvm()` - Matrix-vector multiplication
- Sparse MNA solver for IR drop analysis
- ForwardDiff.jl integration for differentiable simulation


```julia
# Device parameters
params = ThresholdParams(
    R_on = 1e3,      # 1 kΩ LRS
    R_off = 100e3,   # 100 kΩ HRS
    Vth_p = 0.5,     # +0.5V threshold
    Vth_n = -0.5,    # -0.5V threshold
    k_on = 1e4,      # SET rate
    k_off = 1e4      # RESET rate
)

# Simulate with triangular voltage input
result = simulate_memristor(params, (0.0, 2e-3), voltage_func)

# Create crossbar array
xbar = CrossbarArray(32, 32; R_on=1e3, R_off=100e3, R_wire=10.0)

# Perform MVM with IR drop
I_out = simulate_crossbar_mvm_with_ir(xbar, V_input)
```

Stay tuned for the full release!

## References

1. Chua, L. (1971). "Memristor—The missing circuit element." *IEEE Trans. Circuit Theory*
2. Strukov, D. et al. (2008). "The missing memristor found." *Nature*
3. Biolek, Z. et al. (2009). "SPICE model of memristor with nonlinear dopant drift." *Radioengineering*
4. Kvatinsky, S. et al. (2015). "VTEAM: A general model for voltage-controlled memristors." *IEEE TCAS-II*

For more you can check out, the *bibliography.bib* file.


---

Daris Idirene  
Université Paris-Saclay  
Faculty of Sciences - E3A Program