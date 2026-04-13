"""
Connectivity-tax benchmark — simulation utilities.

Compute <M_z> from circuits via:
- Statevector (exact, no shot noise — for circuits without mid-circuit measurement)
- Aer shot simulation (for dynamic circuits with mid-circuit measurement)
"""

from __future__ import annotations

import numpy as np
from qiskit.circuit import QuantumCircuit
from qiskit.quantum_info import Statevector, SparsePauliOp


def magnetization_op(n: int) -> SparsePauliOp:
    """M_z = (1/N) sum_i Z_i."""
    terms = [("Z", [i], 1.0 / n) for i in range(n)]
    return SparsePauliOp.from_sparse_list(terms, num_qubits=n)


# ---------------------------------------------------------------------------
# Layout-aware helpers
# ---------------------------------------------------------------------------

def get_data_qubit_indices(
    circuit: QuantumCircuit,
    n_data: int,
    data_reg_name: str = "d",
) -> list[int]:
    """Determine which physical qubit indices in *circuit* correspond to
    the logical data qubits.

    Strategy (tried in order):
    1. If the circuit carries a transpiler layout, use it to map the first
       *n_data* logical qubits (or the register named *data_reg_name*) to
       their physical positions.
    2. Fall back to ``list(range(n_data))``.

    Parameters
    ----------
    circuit : QuantumCircuit
        Transpiled (or untranspiled) circuit.
    n_data : int
        Number of logical data qubits.
    data_reg_name : str
        Name of the data qubit register in the *original* logical circuit
        (default ``"d"``).  Used to resolve register-based layouts.

    Returns
    -------
    list[int] : physical qubit indices for data qubits, length *n_data*.
    """
    # --- Strategy 0: explicit metadata from build_dynamic_*_physical ---
    meta = getattr(circuit, "metadata", None) or {}
    if "data_physical_qubits" in meta:
        phys = meta["data_physical_qubits"]
        if len(phys) >= n_data:
            return list(phys[:n_data])

    layout = getattr(circuit, "layout", None)
    if layout is None:
        return list(range(n_data))

    # --- Try the initial_layout (TranspileLayout) ---
    try:
        initial = layout.initial_layout
        # initial_layout is a Layout: virtual qubit -> physical qubit
        # We want the physical indices of the first n_data virtual qubits
        # (or the virtual qubits from the data register).

        # Attempt register-based lookup first
        input_circuit = layout.input_qubit_mapping
        # input_qubit_mapping: dict[Qubit -> int] mapping original virtual
        # qubits to their index in the original circuit.
        # initial_layout: Layout mapping virtual qubits of the *transpiled*
        # circuit to physical qubits.

        # Build reverse map: original logical index -> transpiled virtual qubit
        # The TranspileLayout stores:
        #   initial_layout : physical -> virtual (of transpiled circuit)
        #   input_qubit_mapping: original_virtual -> index
        # We need: original logical data qubit index -> physical qubit index

        # Get the original virtual qubits ordered by their logical index
        if input_circuit is not None:
            # input_qubit_mapping: {original_qubit: logical_index}
            orig_qubits_by_idx = sorted(input_circuit.items(), key=lambda x: x[1])

            # Try to identify data qubits by register name
            data_orig_qubits = []
            for qubit, idx in orig_qubits_by_idx:
                reg = _find_register_for_qubit(qubit, data_reg_name)
                if reg is not None:
                    data_orig_qubits.append((qubit, idx))

            # If register lookup didn't work, take the first n_data qubits
            if len(data_orig_qubits) < n_data:
                data_orig_qubits = orig_qubits_by_idx[:n_data]

            # Now map from logical index to physical qubit via the
            # final_layout (if routing changed things) or initial_layout.
            phys_indices = []
            final = getattr(layout, "final_layout", None)
            for qubit, orig_idx in data_orig_qubits[:n_data]:
                # initial_layout maps physical -> virtual (transpiled)
                # We need physical index for this original qubit.
                # The original qubit at index orig_idx corresponds to
                # the virtual qubit at the same index in initial_layout.
                phys = layout.initial_virtual_to_physical(orig_idx)
                # If there's a final_layout, it tracks how routing permuted
                # the qubits.  The net physical qubit is:
                #   initial maps logical -> physical_before_routing
                #   final maps physical_after_routing -> virtual
                # For observable measurement we want the *final* physical
                # position (after all SWAPs).
                if final is not None:
                    # final_layout: virtual -> physical (post-routing)
                    # The qubit that started at phys may have been swapped.
                    # Qiskit's TranspileLayout.routing_permutation() gives
                    # the full permutation.
                    pass  # handled below
                phys_indices.append(phys)

            # Apply routing permutation if available
            if final is not None:
                try:
                    perm = layout.routing_permutation()
                    # perm[i] = j means virtual qubit i ends up at physical j
                    phys_indices = [perm[p] for p in phys_indices]
                except Exception:
                    pass

            if len(phys_indices) == n_data:
                return phys_indices
    except Exception:
        pass

    # --- Simpler fallback: use initial_layout directly ---
    try:
        initial = layout.initial_layout
        # Try to get physical qubit for each logical index
        phys = []
        for i in range(n_data):
            p = initial[i]  # Layout.__getitem__ with int key
            phys.append(p)

        # Apply routing permutation if available
        final = getattr(layout, "final_layout", None)
        if final is not None:
            try:
                perm = layout.routing_permutation()
                phys = [perm[p] for p in phys]
            except Exception:
                pass
        return phys
    except Exception:
        pass

    return list(range(n_data))


def _find_register_for_qubit(qubit, target_reg_name: str):
    """Check if a qubit belongs to a register with the given name."""
    try:
        # Qiskit Qubit objects have a _register attribute
        reg = getattr(qubit, "_register", None)
        if reg is not None and reg.name == target_reg_name:
            return reg
    except Exception:
        pass
    return None


def simulate_mz_statevector(circuit: QuantumCircuit, n_data: int,
                            data_qubits: list[int] | None = None) -> float:
    """Compute exact <M_z> via statevector simulation.

    Works for circuits WITHOUT mid-circuit measurements (ideal, unitary).
    The circuit should have measure=False.

    Parameters
    ----------
    circuit : QuantumCircuit
        Circuit without measurements (use measure=False when building).
    n_data : int
        Number of data qubits.
    data_qubits : list[int] or None
        Physical qubit indices that carry the data.  If None, attempts to
        extract from the circuit's transpiler layout; falls back to
        ``range(n_data)``.

    Returns
    -------
    float : exact <M_z>.
    """
    if data_qubits is None:
        data_qubits = get_data_qubit_indices(circuit, n_data)

    sv = Statevector.from_instruction(circuit)
    n_total = circuit.num_qubits
    terms = [("Z", [q], 1.0 / n_data) for q in data_qubits]
    Mz = SparsePauliOp.from_sparse_list(terms, num_qubits=n_total)
    return sv.expectation_value(Mz).real


def simulate_mz_aer(circuit: QuantumCircuit, n_data: int,
                    shots: int = 50000, noise_model=None,
                    data_qubits: list[int] | None = None) -> float:
    """Compute <M_z> via Aer shot simulation.

    Works for all circuits including dynamic circuits with mid-circuit
    measurements.  The circuit should have measure=True with a 'meas'
    classical register for the final data-qubit measurement.

    Parameters
    ----------
    circuit : QuantumCircuit
        Circuit with final measurements on data qubits.
    n_data : int
        Number of data qubits.
    shots : int
        Number of measurement shots.
    noise_model : optional
        Qiskit Aer NoiseModel.  If None, runs noiseless simulation.
    data_qubits : list[int] or None
        Physical qubit indices that carry the data (used when interpreting
        ``measure_all`` results).  Ignored when the circuit has a 'meas'
        classical register, since those bits already map to data qubits.

    Returns
    -------
    float : estimated <M_z>.
    """
    from qiskit_aer import AerSimulator

    if noise_model is not None:
        sim = AerSimulator(noise_model=noise_model)
    else:
        sim = AerSimulator()

    job = sim.run(circuit, shots=shots)
    counts = job.result().get_counts()

    # Determine whether we have a dedicated 'meas' register
    has_meas_reg = any(creg.name == "meas" for creg in circuit.cregs)

    total = sum(counts.values())
    mz = 0.0
    for bitstring, count in counts.items():
        if has_meas_reg:
            # Dynamic circuits have multi-register keys: "meas_bits mid_bits"
            # 'meas' was added last -> appears first (leftmost)
            parts = bitstring.split()
            meas_str = parts[0] if len(parts) > 1 else bitstring
            # Qiskit convention: rightmost bit = qubit 0
            bits = [int(b) for b in reversed(meas_str)]
            zi_sum = sum(1 - 2 * b for b in bits[:n_data])
        else:
            # measure_all: bitstring covers all physical qubits
            # Qiskit convention: rightmost bit = qubit 0
            all_bits = [int(b) for b in reversed(bitstring)]
            if data_qubits is None:
                data_qubits = get_data_qubit_indices(circuit, n_data)
            zi_sum = sum(1 - 2 * all_bits[q] for q in data_qubits)

        mz += (zi_sum / n_data) * count

    return mz / total


def simulate_mz(circuit: QuantumCircuit, n_data: int,
                shots: int = 50000, noise_model=None,
                method: str = "auto",
                data_qubits: list[int] | None = None) -> float:
    """Compute <M_z> from a circuit, choosing the best method.

    Parameters
    ----------
    circuit : QuantumCircuit
        The circuit to simulate.
    n_data : int
        Number of data qubits.
    shots : int
        Shots for Aer simulation (ignored for statevector).
    noise_model : optional
        Noise model for Aer (forces Aer method).
    method : str
        "auto" (default): use statevector if possible, else Aer.
        "statevector": force statevector (fails if mid-circuit measurements).
        "aer": force Aer shot simulation.
    data_qubits : list[int] or None
        Physical qubit indices that carry the data.  If None, extracted
        automatically from the circuit's transpiler layout (falling back
        to ``range(n_data)``).

    Returns
    -------
    float : <M_z>.
    """
    has_midcircuit_meas = "if_else" in circuit.count_ops()
    has_measure = "measure" in circuit.count_ops()

    if noise_model is not None:
        method = "aer"

    if method == "auto":
        if has_midcircuit_meas:
            method = "aer"
        else:
            method = "statevector"

    if method == "statevector":
        # Need a measurement-free version of the circuit
        if has_measure and not has_midcircuit_meas:
            circ_no_meas = circuit.remove_final_measurements(inplace=False)
            return simulate_mz_statevector(circ_no_meas, n_data,
                                           data_qubits=data_qubits)
        elif not has_measure:
            return simulate_mz_statevector(circuit, n_data,
                                           data_qubits=data_qubits)
        else:
            raise ValueError(
                "Cannot use statevector with mid-circuit measurements. "
                "Use method='aer'."
            )
    else:
        # Aer needs measurements
        if not has_measure:
            circ_meas = circuit.copy()
            circ_meas.measure_all()
            return simulate_mz_aer(circ_meas, n_data, shots, noise_model,
                                   data_qubits=data_qubits)
        return simulate_mz_aer(circuit, n_data, shots, noise_model,
                               data_qubits=data_qubits)
