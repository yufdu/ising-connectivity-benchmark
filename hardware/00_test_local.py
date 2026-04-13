# %% [markdown]
# # Hardware Step 0 — Local Test
#
# Test the full hardware workflow locally using a fake (simulated) backend.
# This verifies that transpilation, submission, and result parsing all
# work before spending real QPU time.
#
# **No IBM account needed for this step.**

# %%
import sys, os
sys.path.insert(0, os.path.abspath(".."))

import numpy as np

from src.model import exact_trotter_magnetization, coupling_pairs
from src.circuits import (
    build_ideal_trotter_circuit,
    build_dynamic_trotter_circuit_physical,
    transpile_dynamic_physical,
)
from src.metrics import extract_metrics

from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager
from qiskit_ibm_runtime.fake_provider import FakeSherbrooke
from qiskit_ibm_runtime import SamplerV2 as Sampler

# %% [markdown]
# ## Parameters

# %%
N_QUBITS = 6
ALPHAS_TEST = [1.0, 4.0]       # small subset for quick test
CUTOFF_FRAC = 0.01
THETA = np.pi / 4
N_STEPS = 1
SHOTS = 1000                   # low shots for speed
SEED = 42
GEOMETRY = "chain"

# %% [markdown]
# ## 1. Set up fake backend
#
# FakeSherbrooke mimics a 127-qubit Eagle processor with realistic
# noise.  The circuits run locally but through the same transpilation
# and execution path as real hardware.

# %%
backend = FakeSherbrooke()
cmap = backend.coupling_map
print(f"Fake backend: {backend.name}")
print(f"  qubits: {backend.num_qubits}")

pm_unitary = generate_preset_pass_manager(
    optimization_level=1, backend=backend, seed_transpiler=SEED,
)

# %% [markdown]
# ## Helper: parse counts from the 'meas' register only

# %%
def get_meas_counts(pub_result):
    """Extract counts from the 'meas' classical register.

    SamplerV2 stores each classical register as a separate attribute
    on pub_result.data.  Dynamic circuits have both 'mid' (ancilla
    measurements) and 'meas' (final data measurements).  We must
    read 'meas' explicitly — never fall back to an arbitrary register,
    as that could silently pick 'mid' and produce garbage Mz values.
    """
    if hasattr(pub_result.data, "meas"):
        return pub_result.data.meas.get_counts()
    if hasattr(pub_result.data, "c"):
        # measure_all() in older Qiskit versions names the register 'c'
        return pub_result.data.c.get_counts()
    available = [a for a in dir(pub_result.data) if not a.startswith("_")]
    raise RuntimeError(
        f"Cannot find 'meas' or 'c' register in result data. "
        f"Available registers: {available}.  Check that circuits "
        f"were built with measure=True."
    )


def mz_from_counts(counts, n_data):
    """Compute <M_z> from measurement counts on data qubits."""
    total = sum(counts.values())
    mz = 0.0
    for bitstring, count in counts.items():
        # Qiskit bitstrings: rightmost bit = qubit 0
        bits = [int(b) for b in reversed(bitstring)]
        mz += sum(1 - 2 * b for b in bits[:n_data]) / n_data * count
    return mz / total

# %% [markdown]
# ## 2. Transpile and inspect circuits

# %%
circuits = []
labels = []

for alpha in ALPHAS_TEST:
    # --- Unitary ---
    circ = build_ideal_trotter_circuit(
        N_QUBITS, alpha, THETA, N_STEPS, measure=True, cutoff_frac=CUTOFF_FRAC, geometry=GEOMETRY, 
    )
    circ_t = pm_unitary.run(circ)
    circuits.append(circ_t)
    labels.append(("unitary", alpha))
    m = extract_metrics(circ_t, "unitary", N_QUBITS, alpha, N_STEPS, THETA)
    print(f"unitary α={alpha}: 2Q_depth={m.two_q_depth}, "
          f"2Q_count={m.two_q_count}, SWAPs={m.swap_count}")

    # --- Dynamic (physical) ---
    circ_d = build_dynamic_trotter_circuit_physical(
        N_QUBITS, alpha, THETA, coupling_map=cmap,
        n_steps=N_STEPS, measure=True, cutoff_frac=CUTOFF_FRAC, geometry=GEOMETRY, 
    )
    circ_dt = transpile_dynamic_physical(circ_d, backend=backend, seed=SEED)
    circuits.append(circ_dt)
    labels.append(("dynamic", alpha))
    md = extract_metrics(circ_dt, "dynamic", N_QUBITS, alpha, N_STEPS, THETA)
    print(f"dynamic α={alpha}: 2Q_depth={md.two_q_depth}, "
          f"2Q_count={md.two_q_count}, SWAPs={md.swap_count}")

print(f"\nTotal circuits: {len(circuits)}")

# %% [markdown]
# ## 3. Run on fake backend (local, with noise)

# %%
sampler = Sampler(mode=backend)

print(f"{'variant':>10} {'alpha':>5} | {'Mz_fake':>8} {'Mz_exact':>8} {'error':>7}")
print("-" * 48)

for circ, (variant, alpha) in zip(circuits, labels):
    job = sampler.run([circ], shots=SHOTS)
    result = job.result()
    pub_result = result[0]

    counts = get_meas_counts(pub_result)
    mz = mz_from_counts(counts, N_QUBITS)
    mz_exact = exact_trotter_magnetization(
        N_QUBITS, alpha, 1.0, THETA, N_STEPS, geometry=GEOMETRY, cutoff_frac=1e-14,
    )
    error = abs(mz - mz_exact)
    print(f"{variant:>10} {alpha:>5.1f} | {mz:>+8.4f} {mz_exact:>+8.4f} {error:>7.4f}")

# %% [markdown]
# ## 4. Check
#
# If this ran without errors, the full hardware workflow is working.
# Expect noisy results (errors of 0.05–0.2) since the fake backend
# includes realistic noise.  On real hardware, results will be similar.
#
# You're ready to proceed to `01a_submit_unitary.py`.
