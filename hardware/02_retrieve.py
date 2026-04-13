# %% [markdown]
# # Hardware Step 2 — Retrieve Results
#
# Load job IDs from the selected submission file, check their status,
# and download results when ready.  Can be run repeatedly until all
# jobs complete.
#
# Works with both `jobs_unitary_*.json` and `jobs_dynamic_*.json` files
# produced by the separate submission scripts.
#
# Saves raw counts to `../results/hw_data/`.

# %%
import sys, os, json, glob, time
sys.path.insert(0, os.path.abspath(".."))

from qiskit_ibm_runtime import QiskitRuntimeService

# %% [markdown]
# ## 1. Find job files and let the user choose

# %%
JOBS_DIR = os.path.join(os.path.dirname(__file__) or ".", "jobs")
HW_DATA_DIR = os.path.join(os.path.dirname(__file__) or ".", "..", "results", "hw_data")
os.makedirs(HW_DATA_DIR, exist_ok=True)

jobs_files = sorted(glob.glob(os.path.join(JOBS_DIR, "jobs_*.json")))
if not jobs_files:
    print("No job files found in hardware/jobs/.")
    print("Run 01a_submit_unitary.py or 01b_submit_dynamic.py first.")
    sys.exit(1)

print("Available job files:")
for i, f in enumerate(jobs_files):
    with open(f) as fh:
        meta = json.load(fh)
    circuit_type = meta.get("circuit_type", "mixed")
    n_jobs = len(meta["jobs"])
    print(f"  [{i}] {os.path.basename(f)}  "
          f"({circuit_type}, {n_jobs} jobs, backend={meta['backend']})")

# To retrieve a specific file, change this index:
#FILE_INDEX = -1  # -1 = most recent
FILE_INDEX = 0 

jobs_file = jobs_files[FILE_INDEX]
print(f"\nRetrieving: {os.path.basename(jobs_file)}")

with open(jobs_file) as f:
    submission = json.load(f)

circuit_type = submission.get("circuit_type", "mixed")
print(f"  Type: {circuit_type}")
print(f"  Backend: {submission['backend']}")
print(f"  Submitted: {submission['submitted_at']}")
print(f"  Jobs: {len(submission['jobs'])}")

# %% [markdown]
# ## 2. Connect and check status

# %%
service = QiskitRuntimeService()

job_records = submission["jobs"]
all_done = True

print(f"\n{'variant':>10} {'alpha':>5} {'job_id':>24} {'status':>12}")
print("-" * 58)

for record in job_records:
    job = service.job(record["job_id"])
    status = job.status()
    record["status"] = status
    if status != "DONE":
        all_done = False
    print(f"{record['variant']:>10} {record['alpha']:>5.1f} "
          f"{record['job_id']:>24} {status:>12}")

if all_done:
    print("\nAll jobs complete!")
else:
    print(f"\nSome jobs still running. Re-run this script later.")
    print(f"Or wait here — polling every 30 seconds...")

# %% [markdown]
# ## 3. Wait for completion (optional)

# %%
if not all_done:
    while True:
        time.sleep(30)
        pending = 0
        for record in job_records:
            if record.get("status") != "DONE":
                job = service.job(record["job_id"])
                status = job.status().name
                record["status"] = status
                if status != "DONE":
                    pending += 1
                else:
                    print(f"  Completed: {record['variant']} α={record['alpha']}")
        if pending == 0:
            print("\nAll jobs complete!")
            break
        print(f"  {pending} jobs still pending...")

# %% [markdown]
# ## Helper: safe counts extraction

# %%
def get_meas_counts(pub_result):
    """Extract counts from the 'meas' classical register only.

    SamplerV2 stores each classical register as a separate attribute.
    Dynamic circuits have both 'mid' and 'meas' — we must read 'meas'
    explicitly.  Never fall back to an arbitrary register.
    """
    if hasattr(pub_result.data, "meas"):
        return pub_result.data.meas.get_counts()
    if hasattr(pub_result.data, "c"):
        return pub_result.data.c.get_counts()
    available = [a for a in dir(pub_result.data) if not a.startswith("_")]
    raise RuntimeError(
        f"Cannot find 'meas' or 'c' register in result data. "
        f"Available registers: {available}"
    )

# %% [markdown]
# ## 4. Download results

# %%
N_QUBITS = submission["parameters"]["n_qubits"]
results = []

for record in job_records:
    job = service.job(record["job_id"])
    result = job.result()
    pub_result = result[0]

    counts = get_meas_counts(pub_result)

    # Compute M_z from the 'meas' register counts
    total = sum(counts.values())
    mz = 0.0
    for bitstring, count in counts.items():
        bits = [int(b) for b in reversed(bitstring)]
        mz += sum(1 - 2 * b for b in bits[:N_QUBITS]) / N_QUBITS * count
    mz /= total

    result_record = {
        "job_id": record["job_id"],
        "variant": record["variant"],
        "alpha": record["alpha"],
        "n_qubits": N_QUBITS,
        "n_steps": record["n_steps"],
        "theta": record["theta"],
        "cutoff_frac": record["cutoff_frac"],
        "shots": record["shots"],
        "backend": record.get("backend", submission["backend"]),
        "mz_hardware": round(mz, 6),
        "counts": counts,
    }
    results.append(result_record)
    print(f"  {record['variant']:>10} α={record['alpha']:.1f}: Mz = {mz:+.4f}")

# %% [markdown]
# ## 5. Save results

# %%
timestamp = submission["submitted_at"]
output_file = os.path.join(
    HW_DATA_DIR, f"hw_results_{circuit_type}_{timestamp}.json"
)

with open(output_file, "w") as f:
    json.dump(results, f, indent=2)

print(f"\nResults saved to: {output_file}")
print(f"Run 03_analyze.py to compare with exact simulation.")
