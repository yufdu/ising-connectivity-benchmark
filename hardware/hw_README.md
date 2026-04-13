# Hardware Submission Guide

Run your connectivity-tax benchmark circuits on real IBM quantum hardware.

## Prerequisites

1. **IBM Quantum account** (free): Sign up at [quantum.cloud.ibm.com](https://quantum.cloud.ibm.com/)

2. **Save your credentials** (one-time setup):
   ```python
   from qiskit_ibm_runtime import QiskitRuntimeService
   QiskitRuntimeService.save_account(
       token="<your-44-character-API-key>",   # from IBM Quantum dashboard
       instance="<your-CRN>",                 # from IBM Quantum > Instances page
       overwrite=True,
   )
   ```
   Your API key is on the [dashboard](https://quantum.cloud.ibm.com/) (click "Copy" under API token). Your CRN (Cloud Resource Name) is on the [Instances](https://quantum.cloud.ibm.com/instances) page — hover over the instance and click to copy.

3. **Install dependencies**:
   ```bash
   pip install qiskit-ibm-runtime
   ```

4. **Run notebooks 01–03 first** to generate the local simulation results.

## Workflow

### Step 0: Local test (no account needed)
```bash
python 00_test_local.py
```
Verifies that transpilation, execution, and result parsing all work for both unitary and dynamic circuits. 

### Step 1a: Submit unitary circuits
```bash
python 01a_submit_unitary.py
```
Transpiles ideal Trotter circuits with SWAP routing and submits them. Saves job IDs to `jobs/jobs_unitary_YYYYMMDD_HHMMSS.json`. 

### Step 1b: Submit dynamic circuits
```bash
python 01b_submit_dynamic.py
```
Builds dynamic circuits directly on the backend's physical topology and submits them. Saves job IDs to `jobs/jobs_dynamic_YYYYMMDD_HHMMSS.json`.

The unitary and dynamic submissions are independent. You can run 1a now, iterate on the dynamic circuit design, and run 1b later without affecting the unitary results.

### Step 2: Retrieve results
```bash
python 02_retrieve.py
```
Lists all job files in `jobs/` and retrieves the most recent one. Change `FILE_INDEX` in the script to retrieve a different submission. Can be run repeatedly — safe to re-run.

Saves results to `../results/hw_data/hw_results_{type}_{timestamp}.json`.

### Step 3: Analyze
```bash
python 03_analyze.ipynb
```
Loads **all** result files in `results/hw_data/`, merges them, and produces comparison plots. Works whether you have just unitary results, just dynamic, or both. Runs entirely offline — no IBM account needed.

Saves plot to `../results/hw_data/hardware_results.png`.

## File structure

```
hardware/
├── 00_test_local.py          # Step 0: local test with fake backend
├── 01a_submit_unitary.py     # Step 1a: submit unitary (SWAP-routed)
├── 01b_submit_dynamic.py     # Step 1b: submit dynamic (ancilla-assisted)
├── 02_retrieve.py            # Step 2: poll status + download results
├── 03_analyze.ipynb          # Step 3: compare with exact, make plots
├── jobs/                     # auto-created: job ID records
│   ├── jobs_unitary_xxxxxxxx_xxxxxx.json
│   └── jobs_dynamic_xxxxxxxx_xxxxxx.json
└── README.md                 # this file

results/hw_data/              # auto-created: hardware measurement data
    ├── hw_results_unitary_xxxxxxxx_xxxxxx.json
    ├── hw_results_dynamic_xxxxxxxx_xxxxxx.json
    └── hardware_results.png
```

## Circuit parameters and evolution time

The hardware scripts use the following default settings:

| Parameter     | Value                     | Description                              |
| ------------- | ------------------------- | ---------------------------------------- |
| `N_QUBITS`    | 6                         | Data qubits (1D chain geometry)          |
| `N_STEPS`     | 1                         | Number of first-order Trotter steps      |
| `THETA`       | $\pi/4$                   | Kick angle per Trotter step              |
| `CUTOFF_FRAC` | 0.01                      | Keep coupling pairs with J ≥ 1% of J_max |
| `SHOTS`       | 4000                      | Measurement shots per circuit            |
| `ALPHAS`      | [0.5, 1.0, 2.0, 4.0, 8.0] | Power-law exponents scanned              |

In this code, each Trotter step applies an X-kick $R_x(\theta)$ followed by ZZ interactions $R_{zz}(J_{ij})$. The kick relates to the continuous-time transverse field $h$ and time step $dt$ via $R_x(\theta) = e^{−i \theta X / 2} = e^{−i h dt X}$, so $\theta = 2 h dt$. With the defaults $\theta = \pi/4$ and $h = 1.0$, each step corresponds to $dt = \pi/8 ≈ 0.393$ in natural units. For `N_STEPS = 1`, the total evolution time is $T = dt = \pi/8$.

To increase the evolution time (or number of Trotter steps), edit `N_STEPS` in the submit scripts. The simulation notebooks (01–04) sweep up to `MAX_STEPS = 8` to study how the connectivity tax scales with depth. Notebook 04 shows phase diagram where the 2Q depth exceeds typical coherence budgets.

## Job ID files

Each submission creates a JSON file in `jobs/` containing:
- `submitted_at`: timestamp
- `backend`: which IBM processor was used
- `circuit_type`: `"unitary"` or `"dynamic"`
- `parameters`: N, α values, θ, cutoff, shots
- `jobs`: list of `{job_id, variant, alpha, ...}` records

These files are your receipt. If your Python session crashes after submission, you can recover results using the job IDs — either through `02_retrieve.py` or manually:
```python
from qiskit_ibm_runtime import QiskitRuntimeService
service = QiskitRuntimeService()
job = service.job("your-job-id-here")
print(job.status())
result = job.result()
```

You can also check job status on the [IBM Quantum Workloads](https://quantum.cloud.ibm.com/workloads) page.

## IBM Quantum Open Plan limits

The free Open Plan provides **10 minutes of QPU time every 28 days**. Our benchmark (5 unitary + 5 dynamic circuits × 4000 shots each) uses roughly 1–3 minutes of QPU time, well within the budget. Queue wait time (which can be minutes to hours) does not count against your allocation.

## Troubleshooting

**"Backend does not support control-flow"**: Not all backends support `if_else`. Try a specific backend known to support it:
```python
backend = service.backend("ibm_sherbrooke")
```

**"No backends available"**: The free plan may have limited backend availability. Try again later, or specify a backend directly.

**Job stuck in queue**: Normal on the free plan. Check the [IBM Quantum dashboard](https://quantum.cloud.ibm.com/workloads) for estimated wait times. You can close your terminal and retrieve results later.

**"Credentials not found"**: Re-run the `save_account` step in the Prerequisites section above.
