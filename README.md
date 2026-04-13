# Connectivity Tax: Benchmarking Long-Range Ising Simulation on Heavy-Hex Topology

## Motivation

IBM's superconducting processors use a fixed heavy-hex topology, which forces routing overhead for non-local interactions.  IBM's own [dynamic circuit tutorial](https://quantum.cloud.ibm.com/docs/en/tutorials/dc-hex-ising) shows that mid-circuit measurement and feedforward can bypass some of this overhead, but only when the interaction graph matches the hardware ancilla layout. 

This exercise investigates the heavy-hex connectivity cost for simulating a generic long-range spin model, and explores the behavior of dynamic circuit method beyond topology-matched cases. 

## Model

Long-range Ising Hamiltonian (this version: 1D chain with N = 6 spins): $ H = \sum_{i<j} J_{ij} \, Z_i Z_j + h \sum_{i} X_i $, where $J_{ij} = 1/|i−j|^\alpha$  is a power-law coupling, truncated at a tunable cutoff (set to 1% of max).

- Large $\alpha$: short range interaction (only nearest-neighbor pairs survive, trivial on heavy-hex)
- Small $\alpha$: long range interaction (has long-range pairs, high connectivity demand)

Observable: average magnetization $M_z = (1/N) \sum_i \langle Z_i\rangle$

## Three Benchmark Variants

| Variant               | Connectivity | Routing method                                                                                                                                                                     |
| --------------------- | ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Ideal circuit**     | All-to-all   | None (direct Rzz)                                                                                                                                                                  |
| **Heavy-hex unitary** | Heavy-hex    | SWAP insertion by transpiler                                                                                                                                                       |
| **Heavy-hex dynamic** | Heavy-hex    | - NN pairs: ancilla-assisted ZZ (mid-circuit measurement + feedforward, with each ancilla pinned to its hardware bridge qubit)<br>- Long-range pairs: SWAP insertion by transpiler |

## Geometry and qubit number

(Upcoming...)

This code intends to support both `"chain"` (1D) and `"honeycomb"` (hexagonal lattice) geometries, controlled by a single `geometry` parameter. To switch, set `GEOMETRY = "honeycomb"`.

Acceptable number of qubit is determined by the honeycomb structure:
- N_QUBITS=6  -> hex_lattice(1,1): one hexagonal ring
- N_QUBITS=10 -> hex_lattice(1,2) or (2,1): two fused hexagons
- N_QUBITS=16 -> hex_lattice(2,2): four fused hexagons
etc. 

**WARNING: not fully implemented in this version.**

## Repository Structure

```
connectivity-tax/
├── README.md
├── requirements.txt
├── src/
│   ├── __init__.py
│   ├── model.py        # Hamiltonian, couplings, exact diagonalisation
│   ├── circuits.py     # Circuit builders for all 3 variants
│   ├── metrics.py      # Metric extraction and serialisation
│   ├── simulate.py     # Aer / statevector simulation helpers
│   └── analysis.py     # Plotting utilities
├── notebooks/
│   ├── 01_model_and_ideal.ipynb       # Define model, build ideal circuits
│   ├── 02_heavyhex_unitary.ipynb      # SWAP routing variant
│   ├── 03_heavyhex_dynamic.ipynb      # Dynamic circuit variant
│   ├── 04_benchmark_analysis.ipynb    # All plots and phase diagrams
│   └── 05_circuit_verification.ipynb  # Verify circuits: noiseless simulation
│   
├── hardware/                       # IBM hardware submission (separate workflow)
│   ├── hw_README.md                # Setup guide and troubleshooting
│   ├── 00_test_local.py            # Test on fake backend (no account needed)
│   ├── 01a_submit_unitary.py       # Transpile + submit unitary, save job IDs
│   ├── 01b_submit_dynamic.py       # Transpile + submit dynamic, save job IDs
│   ├── 02_retrieve.py              # Poll status + download results
│   ├── 03_analyze.ipynb            # Compare hardware vs exact
│   └── jobs/                       # Auto-created: saved job ID files
└── results/
    ├── ideal_metrics.json
    ├── unitary_metrics.json
    ├── dynamic_metrics.json
    └── hw_data/                    # Auto-created: hardware measurement data
```

## Running

```bash
pip install -r requirements.txt
```

- Generate all metrics (run notebooks 01-03 in order)
- Analysis (run notebook 04, which loads results from 01-03)
- Verify circuits match exact simulation (recommended before hardware)
- Hardware (requires IBM Quantum account — see hardware/README.md)

Notebooks use `# %%` cell markers — open them directly in VS Code or JupyterLab as interactive notebooks.
## IBM Quantum Hardware (Optional)

Running on real hardware requires a free IBM Quantum account:

1. Sign up at [quantum.cloud.ibm.com](https://quantum.cloud.ibm.com/)
2. Copy your API key from the dashboard
3. Copy your instance CRN from the [Instances](https://quantum.cloud.ibm.com/instances) page
4. Save credentials (one-time):
   ```python
   from qiskit_ibm_runtime import QiskitRuntimeService
   QiskitRuntimeService.save_account(token="<API-key>", instance="<CRN>")
   ```

The free Open Plan provides 10 minutes of QPU time every 28 days — enough for our benchmark. See `hardware/README.md` for the full guide, including troubleshooting and how to recover results from job IDs.

## References

- IBM Quantum, ["Simulation of kicked Ising Hamiltonian with dynamic circuits"](https://quantum.cloud.ibm.com/docs/en/tutorials/dc-hex-ising) 