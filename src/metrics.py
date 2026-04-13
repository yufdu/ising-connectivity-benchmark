"""
Connectivity-tax benchmark — metric extraction.

Extract compile-time metrics from transpiled circuits:
- total depth, 2Q depth, 2Q gate count, SWAP count
- estimated circuit duration (if scheduling info available)
"""

from __future__ import annotations

import json
from dataclasses import dataclass, asdict
from pathlib import Path

from qiskit.circuit import QuantumCircuit


@dataclass
class CircuitMetrics:
    """Container for benchmark metrics of a single circuit."""
    variant: str          # "ideal", "unitary", "dynamic"
    n_qubits: int         # logical qubit count
    alpha: float          # power-law exponent
    n_steps: int          # Trotter steps
    theta: float          # kick angle
    geometry: str = "chain"  # "chain" or "honeycomb"

    total_depth: int = 0
    two_q_depth: int = 0
    two_q_count: int = 0
    swap_count: int = 0
    total_gate_count: int = 0
    cx_count: int = 0
    measure_count: int = 0
    reset_count: int = 0

    # Estimated duration in microseconds (from scheduling, if available)
    estimated_duration_us: float | None = None

    # Observable result from noiseless simulation
    mz_exact: float | None = None
    mz_sim: float | None = None

    def to_dict(self) -> dict:
        return asdict(self)


def extract_metrics(
    circuit: QuantumCircuit,
    variant: str,
    n_qubits: int,
    alpha: float,
    n_steps: int,
    theta: float,
    geometry: str = "chain",
) -> CircuitMetrics:
    """Extract compile-time metrics from a (transpiled) circuit.

    Parameters
    ----------
    circuit : QuantumCircuit
        The circuit to analyse (should be transpiled to a basis gate set).
    variant : str
        Label: "ideal", "unitary", or "dynamic".
    n_qubits, alpha, n_steps, theta : problem parameters.

    Returns
    -------
    CircuitMetrics
    """
    m = CircuitMetrics(
        variant=variant,
        n_qubits=n_qubits,
        alpha=alpha,
        n_steps=n_steps,
        theta=theta,
        geometry=geometry,
    )

    m.total_depth = circuit.depth()
    m.total_gate_count = circuit.size()

    # Count by operation name
    ops = circuit.count_ops()
    two_q_gates = {"cx", "cz", "ecr", "rzz", "rxx", "ryy", "cp", "swap", "iswap"}
    m.two_q_count = sum(v for k, v in ops.items() if k.lower() in two_q_gates)
    m.swap_count = ops.get("swap", 0)
    m.cx_count = ops.get("cx", 0) + ops.get("ecr", 0)
    m.measure_count = ops.get("measure", 0)
    m.reset_count = ops.get("reset", 0)

    # Two-qubit depth: depth considering only 2Q gates
    # We build a filtered circuit to compute this
    try:
        from qiskit.converters import circuit_to_dag, dag_to_circuit
        dag = circuit_to_dag(circuit)
        # Count longest path through 2Q gates
        two_q_depth = 0
        for layer in dag.layers():
            has_2q = any(
                len(node.qargs) >= 2 and node.op.name.lower() in two_q_gates
                for node in layer["graph"].op_nodes()
            )
            if has_2q:
                two_q_depth += 1
        m.two_q_depth = two_q_depth
    except Exception:
        m.two_q_depth = -1  # fallback

    # Try to get scheduled duration
    try:
        if circuit.duration is not None:
            # Duration is in dt units; convert assuming dt ~ 0.222 ns (IBM default)
            m.estimated_duration_us = circuit.duration * 0.222e-3
    except Exception:
        pass

    return m


def save_metrics(metrics_list: list[CircuitMetrics], path: str | Path):
    """Save a list of CircuitMetrics to JSON."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = [m.to_dict() for m in metrics_list]
    with open(path, "w") as f:
        json.dump(data, f, indent=2)


def load_metrics(path: str | Path) -> list[CircuitMetrics]:
    """Load CircuitMetrics from JSON."""
    with open(path) as f:
        data = json.load(f)
    return [CircuitMetrics(**d) for d in data]
