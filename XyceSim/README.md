# XyceSim

Julia interface for Xyce simulator via Jyce wrapper for memristor crossbar simulation.

## Features

- Xyce simulator integration through Jyce wrapper
- Crossbar netlist generation for both ADMS plugin and subcircuit approaches
- Simulation execution and result postprocessing
- Validation plotting and analysis tools
- Support for memristor device models compatible with Xyce

## Installation

First, ensure you have Jyce properly installed and configured:
```bash
# Install Jyce (if not already installed)
julia -e 'using Pkg; Pkg.add("Jyce")'

# Configure Jyce with your Xyce installation
# Set environment variables JYCE_XYCESOLVER_JULIA_LIB or JYCE_XYCESOLVER_ROOT
# as needed for your system
```

Then add XyceSim:
```julia
using Pkg
Pkg.add("XyceSim")
```

## Usage

```julia
using XyceSim

# Create simulator instance
sim = XyceSim.create_simulator()

# Generate crossbar netlist using ADMS memristor plugin
params = XyceSim.MemristorParams(Ron=1e3, Roff=100e3, D=3e-9, uv=1e-15)
crossbar = XyceSim.create_crossbar_array(4, 4, params=params)

# Write netlist to file
XyceSim.write_netlist(crossbar, "crossbar_netlist.cir")

# Run simulation
prn_file = XyceSim.run_crossbar(sim, crossbar)

# Generate plots
results = XyceSim.run_and_plot_crossbar(sim, crossbar)
```

## Documentation

See the `docs/` directory for detailed documentation generation scripts.

## Testing

Run tests with:
```julia
using Pkg
Pkg.test("XyceSim")
```

## Dependencies

- Jyce: Julia wrapper for Xyce parallel electronic simulator
- Plots: Visualization
- Statistics: Statistical functions
- LinearAlgebra: Linear algebra operations
- Random: Random number generation
- Printf: Formatted printing