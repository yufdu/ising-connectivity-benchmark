# Connectivity Tax: Benchmarking Heavy-Hex Topology for Long-Range Ising Simulation

## Motivation

Recent work on Shor's algorithm with reconfigurable neutral atoms ([arXiv:2603.28627](https://arxiv.org/abs/2603.28627)) highlights the advantage of all-to-all connectivity for executing non-local quantum operations.  IBM's superconducting processors use a fixed heavy-hex topology, which forces routing overhead for non-local interactions.  IBM's own [dynamic circuit tutorial](https://quantum.cloud.ibm.com/docs/en/tutorials/dc-hex-ising) shows that mid-circuit measurement and feedforward can bypass some of this overhead — but only when the interaction graph matches the hardware ancilla layout.

**This project asks**: How much does heavy-hex connectivity cost for simulating a generic long-range spin model, and does the dynamic circuit trick generalise beyond topology-matched cases?

## Model

Long-range Ising Hamiltonian on a 1D chain of N = 6 spins:

```
H = Σ_{i<j} J_{ij} Z_i Z_j  +  h Σ_i X_i
```

where `J_{ij} = 1/|i−j|^α` (power-law coupling, truncated at 10% of max).

- **α → large**: only nearest-neighbour pairs survive → trivial on heavy-hex
- **α → small**: many long-range pairs → high connectivity demand

Observable: average magnetisation `M_z = (1/N) Σ_i ⟨Z_i⟩`

## Three Benchmark Variants

| Variant | Connectivity | Routing method |
|---------|-------------|----------------|
| **Ideal** | All-to-all | None (direct Rzz) |
| **Heavy-hex unitary** | Heavy-hex (57 qubits) | SWAP insertion by transpiler |
| **Heavy-hex dynamic** | Heavy-hex (57 qubits) | Ancilla + mid-circuit measurement + feedforward |

## Key Results (N=6, 1 Trotter step)

| α | Active pairs | Ideal 2Q depth | Unitary 2Q depth | Dynamic 2Q depth | Uni tax | Dyn tax |
|---|-------------|----------------|-------------------|-------------------|---------|---------|
| 0.5 | 15 | 9 | 16 | 26 | 1.78× | 2.89× |
| 1.5 | 14 | 8 | 18 | 24 | 2.25× | 3.00× |
| 2.0 | 12 | 7 | 12 | 22 | 1.71× | 3.14× |
| 3.0 | 9 | 6 | 8 | 12 | 1.33× | 2.00× |
| 4.0 | 5 | 5 | 5 | 8 | 1.00× | 1.60× |

## Key Findings

1. **The connectivity tax is real and measurable.** SWAP routing on heavy-hex inflates 2Q depth by up to 2.25× for strongly nonlocal interactions.

2. **Dynamic circuits do NOT recover the overhead for arbitrary long-range chains.** They are consistently 1.3–1.8× more expensive than SWAP routing for our model.

3. **Why the IBM tutorial works but this doesn't:** The tutorial's hexagonal Ising model maps perfectly onto the heavy-hex bridge structure — every pair's ancilla is already physically adjacent to both data qubits.  Our 1D long-range chain lacks this topology match, so ancillas need their own routing.

4. **The feasibility boundary** (2Q depth ≤ 50) limits useful simulation to α ≳ 2.0 at 8 Trotter steps.  Below this, reconfigurable all-to-all connectivity is needed.

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
│   ├── 01_model_and_ideal.py       # Define model, build ideal circuits
│   ├── 02_heavyhex_unitary.py      # SWAP routing variant
│   ├── 02a_verification.py        # Verify all circuits via noiseless simulation
│   ├── 03_heavyhex_dynamic.py      # Dynamic circuit variant
│   ├── 04_benchmark_analysis.py    # All plots and phase diagrams
│   └── 06_crossover_analysis.py    # When do dynamic circuits win?
├── hardware/                       # IBM hardware submission (separate workflow)
│   ├── README.md                   # Setup guide and troubleshooting
│   ├── test_local.py               # Test on fake backend (no account needed)
│   ├── submit.py                   # Transpile + submit, save job IDs
│   ├── retrieve.py                 # Poll status + download results
│   ├── analyze.py                  # Compare hardware vs exact
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

# Generate all metrics (notebooks 01-03 can run in any order)
cd notebooks
python 01_model_and_ideal.py
python 02_heavyhex_unitary.py
python 03_heavyhex_dynamic.py

# Verify circuits match exact simulation (recommended before hardware)
python 02a_verification.py

# Analysis (loads results from 01-03)
python 04_benchmark_analysis.py

# Crossover analysis (transpilation only, no simulation)
python 06_crossover_analysis.py

# Hardware (requires IBM Quantum account — see hardware/README.md)
cd ../hardware
python test_local.py        # test locally first
python submit.py            # submit to IBM
python retrieve.py          # download results when ready
python analyze.py           # compare with exact
```

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

- M. Cain et al., "Shor's algorithm is possible with as few as 10,000 reconfigurable atomic qubits," [arXiv:2603.28627](https://arxiv.org/abs/2603.28627) (2026).
- IBM Quantum, ["Simulation of kicked Ising Hamiltonian with dynamic circuits"](https://quantum.cloud.ibm.com/docs/en/tutorials/dc-hex-ising) — the interaction-topology-matched case where dynamic circuits excel.
- E. Bäumer et al., "Efficient Long-Range Entanglement using Dynamic Circuits," [arXiv:2308.13065](https://arxiv.org/abs/2308.13065) (2023).
- IBM Quantum, ["Run your first circuit on hardware"](https://quantum.cloud.ibm.com/docs/en/guides/hello-world) — hello-world tutorial for IBM Quantum.
