# %% [markdown]
# # Hardware Step 1a — Submit Unitary Jobs
#
# Transpile ideal Trotter circuits (SWAP-routed) for a real IBM backend
# and submit them.  Saves job IDs to `jobs/` for later retrieval.
#
# **Requires IBM Quantum account.  Run 00_test_local.py first.**

# %%
import sys, os, json
from datetime import datetime

sys.path.insert(0, os.path.abspath(".."))

import numpy as np

from src.model import coupling_pairs
from src.circuits import build_ideal_trotter_circuit
from src.metrics import extract_metrics

from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler

# %% [markdown]
# ## Parameters

# %%
N_QUBITS = 6
ALPHAS_HW = [0.5, 1.0, 2.0, 4.0, 8.0]
CUTOFF_FRAC = 0.01
THETA = np.pi / 4
N_STEPS = 1
SHOTS = 4000
SEED = 42
GEOMETRY = "chain"

JOBS_DIR = os.path.join(os.path.dirname(__file__) or ".", "jobs")
os.makedirs(JOBS_DIR, exist_ok=True)

# %% [markdown]
# ## 1. Connect and select backend

# %%
# First-time setup — save your credentials (uncomment and fill in):
# QiskitRuntimeService.save_account(
#     token="<your-44-character-API-key>",   # from IBM Quantum dashboard
#     instance="<your-CRN>",                 # from IBM Quantum > Instances
#     overwrite=True,
# )

service = QiskitRuntimeService()
backend = service.least_busy(simulator=False, operational=True, min_num_qubits=27)
print(f"Backend: {backend.name}  ({backend.num_qubits} qubits)")

# %% [markdown]
# ## 2. Transpile circuits

# %%
pm = generate_preset_pass_manager(
    optimization_level=1, backend=backend, seed_transpiler=SEED,
)

circuits = []
labels = []

for alpha in ALPHAS_HW:
    circ = build_ideal_trotter_circuit(
        N_QUBITS, alpha, THETA, N_STEPS, measure=True, cutoff_frac=CUTOFF_FRAC, geometry=GEOMETRY
    )
    circ_t = pm.run(circ)
    circuits.append(circ_t)
    labels.append({"variant": "unitary", "alpha": alpha})

    m = extract_metrics(circ_t, "unitary", N_QUBITS, alpha, N_STEPS, THETA)
    print(f"unitary α={alpha}: 2Q_depth={m.two_q_depth}, "
          f"2Q_count={m.two_q_count}, SWAPs={m.swap_count}")

print(f"\nTotal circuits to submit: {len(circuits)}")

# %% [markdown]
# ## 3. Submit jobs

# %%
sampler = Sampler(mode=backend)

job_records = []
for circ, label in zip(circuits, labels):
    job = sampler.run([circ], shots=SHOTS)
    record = {
        "job_id": job.job_id(),
        "variant": label["variant"],
        "alpha": label["alpha"],
        "n_qubits": N_QUBITS,
        "n_steps": N_STEPS,
        "theta": THETA,
        "cutoff_frac": CUTOFF_FRAC,
        "geometry": GEOMETRY, 
        "shots": SHOTS,
        "backend": backend.name,
    }
    job_records.append(record)
    print(f"  Submitted unitary α={label['alpha']}: {job.job_id()}")

# %% [markdown]
# ## 4. Save job IDs

# %%
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
jobs_file = os.path.join(JOBS_DIR, f"jobs_unitary_{timestamp}.json")

submission_record = {
    "submitted_at": timestamp,
    "backend": backend.name,
    "circuit_type": "unitary",
    "parameters": {
        "n_qubits": N_QUBITS,
        "alphas": ALPHAS_HW,
        "cutoff_frac": CUTOFF_FRAC,
        "geometry": GEOMETRY, 
        "theta": THETA,
        "n_steps": N_STEPS,
        "shots": SHOTS,
    },
    "jobs": job_records,
}

with open(jobs_file, "w") as f:
    json.dump(submission_record, f, indent=2)

print(f"\nJob IDs saved to: {jobs_file}")
print(f"Run 02_retrieve.py to check status and download results.")
