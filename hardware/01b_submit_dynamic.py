# %% [markdown]
# # Hardware Step 1b — Submit Dynamic Jobs
#
# Build dynamic Trotter circuits directly on the backend's physical
# topology, transpile, and submit.  This script is separate from the
# unitary submission so you can iterate on the dynamic circuits without
# re-running the unitary baseline.
#
# **Requires IBM Quantum account.  Run 00_test_local.py first.**
#
# **Important**: Dynamic circuits use control flow (`if_else`).  The
# backend must support this, and fractional gates must be disabled
# (they are incompatible with control-flow instructions).

# %%
import sys, os, json
from datetime import datetime

sys.path.insert(0, os.path.abspath(".."))

import numpy as np

from src.model import coupling_pairs
from src.circuits import (
    build_dynamic_trotter_circuit_physical,
    transpile_dynamic_physical,
)
from src.metrics import extract_metrics

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
#
# We validate that the chosen backend supports dynamic circuits
# (control flow).  If not, the script fails immediately rather than
# silently submitting an incomplete experiment.

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

# Validate dynamic circuit support
target = backend.target
supported_ops = target.operation_names
if "if_else" not in supported_ops:
    raise RuntimeError(
        f"Backend {backend.name} does not support control-flow (if_else). "
        f"Dynamic circuits require a backend with control-flow support. "
        f"Try selecting a specific backend known to support it, e.g.:\n"
        f"  backend = service.backend('ibm_sherbrooke')"
    )
print(f"  Control-flow support: ✓ (if_else in target)")

# %% [markdown]
# ## 2. Build and transpile circuits
#
# The dynamic circuits are built directly on the backend's physical
# coupling map.  NN ancilla blocks are placed on hardware bridge
# qubits — no SWAPs needed for those.  Only long-range Rzz gates
# require routing.
#
# Any transpilation failure is **fatal** — we do not silently skip
# circuits, as that would produce an incomplete benchmark.

# %%
cmap = backend.coupling_map

circuits = []
labels = []

for alpha in ALPHAS_HW:
    circ_d = build_dynamic_trotter_circuit_physical(
        N_QUBITS, alpha, THETA, coupling_map=cmap,
        n_steps=N_STEPS, measure=True, cutoff_frac=CUTOFF_FRAC, geometry=GEOMETRY
    )

    # Transpile using the backend's full target (basis gates, timing,
    # calibrations) so the output is ISA-compliant.  Identity layout
    # ensures NN ancilla blocks stay on their bridge qubits.
    circ_dt = transpile_dynamic_physical(circ_d, backend=backend, seed=SEED)

    circuits.append(circ_dt)
    labels.append({"variant": "dynamic", "alpha": alpha})

    md = extract_metrics(circ_dt, "dynamic", N_QUBITS, alpha, N_STEPS, THETA)
    print(f"dynamic α={alpha}: 2Q_depth={md.two_q_depth}, "
          f"2Q_count={md.two_q_count}, SWAPs={md.swap_count}")

print(f"\nAll {len(circuits)} circuits transpiled successfully.")

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
    print(f"  Submitted dynamic α={label['alpha']}: {job.job_id()}")

# %% [markdown]
# ## 4. Save job IDs

# %%
timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
jobs_file = os.path.join(JOBS_DIR, f"jobs_dynamic_{timestamp}.json")

submission_record = {
    "submitted_at": timestamp,
    "backend": backend.name,
    "circuit_type": "dynamic",
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
