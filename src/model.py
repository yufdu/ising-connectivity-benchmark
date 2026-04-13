"""
Connectivity-tax benchmark — model definitions.

Long-range Ising Hamiltonian:
    H = sum_{i<j} J_{ij} Z_i Z_j  +  h sum_i X_i
with power-law coupling  J_{ij} = 1 / d(i,j)^alpha,
where d(i,j) is the distance on the underlying geometry (1D chain or honeycomb).

Observable:  M_z = (1/N) sum_i <Z_i>

Relation to a continuous-time transverse field h over one fixed step dt:
    Rx(THETA) = exp(-i THETA X / 2) = exp(-i h dt X)
so
    THETA = 2 h dt.

Since this project keeps dt fixed and scans THETA directly, THETA is the
primary control parameter.

"""

import numpy as np
from scipy.linalg import expm
from qiskit.quantum_info import SparsePauliOp, Statevector


# ---------------------------------------------------------------------------
# Geometry: distance matrices
# ---------------------------------------------------------------------------

def _chain_distances(n: int) -> np.ndarray:
    """Distance matrix for a 1D chain: d(i,j) = |i - j|."""
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            D[i, j] = abs(i - j)
    return D


def _honeycomb_distances(n: int = 6) -> np.ndarray:
    """Distance matrix for a honeycomb (hexagonal) lattice.

    Uses Qiskit's from_hexagonal_lattice to build the graph, then computes
    all-pairs shortest paths.

    Supported sizes:
        n=6  -> hex_lattice(1,1): one hexagonal ring
        n=10 -> hex_lattice(1,2) or (2,1): two fused hexagons
        n=16 -> hex_lattice(2,2): four fused hexagons

    For n=6, the graph is:  0-1-2
                             |   |
                             3-4-5
    with edges: (0,1),(1,2),(0,3),(3,4),(4,5),(2,5)
    """
    from qiskit.transpiler import CouplingMap
    import rustworkx as rx

    # Map n to (rows, cols)
    size_map = {6: (1, 1), 10: (1, 2), 16: (2, 2), 12: (2, 1)}
    if n not in size_map:
        # Try to find a configuration that gives n qubits
        for r in range(1, 6):
            for c in range(1, 6):
                cmap = CouplingMap.from_hexagonal_lattice(r, c, bidirectional=False)
                if cmap.size() == n:
                    size_map[n] = (r, c)
                    break
            if n in size_map:
                break
    if n not in size_map:
        raise ValueError(
            f"Cannot build a honeycomb lattice with exactly {n} qubits. "
            f"Supported sizes: 6, 10, 12, 16, or check from_hexagonal_lattice."
        )

    rows, cols = size_map[n]
    cmap = CouplingMap.from_hexagonal_lattice(rows, cols, bidirectional=False)
    g = cmap.graph.to_undirected(multigraph=False)

    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        lengths = rx.dijkstra_shortest_path_lengths(g, i, lambda _: 1.0)
        for j, d in lengths.items():
            D[i][j] = d
    return D


def distance_matrix(n: int, geometry: str = "chain") -> np.ndarray:
    """Return the distance matrix for the given geometry.

    Parameters
    ----------
    n : int               Number of qubits.
    geometry : str        "chain" or "honeycomb".
    """
    if geometry == "chain":
        return _chain_distances(n)
    elif geometry == "honeycomb":
        return _honeycomb_distances(n)
    else:
        raise ValueError(f"Unknown geometry: {geometry!r}. Use 'chain' or 'honeycomb'.")


# ---------------------------------------------------------------------------
# Coupling matrix
# ---------------------------------------------------------------------------

def coupling_matrix(n: int, alpha: float, geometry: str = "chain") -> np.ndarray:
    """Return the N×N coupling matrix  J[i,j] = 1/d(i,j)^alpha  (upper-triangular).

    Parameters
    ----------
    n : int
        Number of qubits (spins).
    alpha : float
        Power-law exponent.  alpha -> inf  gives nearest-neighbour only.
    geometry : str
        "chain" (1D) or "honeycomb" (hexagonal lattice).

    Returns
    -------
    J : np.ndarray of shape (n, n)
        Upper-triangular coupling matrix.  J[i,j] > 0 for i < j.
    """
    D = distance_matrix(n, geometry)
    J = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            d = D[i, j]
            if d > 0:
                J[i, j] = 1.0 / d ** alpha
    return J


def coupling_pairs(n: int, alpha: float, cutoff_frac: float = 0.10,
                   geometry: str = "chain"):
    """Return list of (i, j, J_ij) for pairs with J_ij >= cutoff_frac * J_max.

    The cutoff makes alpha physically meaningful for circuit structure:
    large alpha -> weak long-range couplings fall below cutoff -> fewer pairs.

    Parameters
    ----------
    n : int           Number of spins.
    alpha : float     Power-law exponent.
    cutoff_frac : float
        Keep pairs with J_ij >= cutoff_frac * max(J).  Default 0.10 (10%).
    geometry : str    "chain" or "honeycomb".

    Returns
    -------
    List of (i, j, J_ij) sorted by descending coupling strength.
    """
    J = coupling_matrix(n, alpha, geometry=geometry)
    J_max = J.max()
    threshold = cutoff_frac * J_max if J_max > 0 else 0.0
    pairs = []
    for i in range(n):
        for j in range(i + 1, n):
            if J[i, j] >= threshold:
                pairs.append((i, j, J[i, j]))
    pairs.sort(key=lambda x: -x[2])
    return pairs


# ---------------------------------------------------------------------------
# Hamiltonian as SparsePauliOp
# ---------------------------------------------------------------------------

def build_hamiltonian(n: int, alpha: float, h: float = 1.0,
                     geometry: str = "chain") -> SparsePauliOp:
    """Build H = sum_{i<j} J_{ij} Z_i Z_j  +  h sum_i X_i.

    Parameters
    ----------
    n : int
        Number of qubits.
    alpha : float
        Power-law exponent for couplings.
    h : float
        Transverse field strength.

    Returns
    -------
    SparsePauliOp
    """
    terms = []
    # ZZ terms
    J = coupling_matrix(n, alpha, geometry=geometry)
    for i in range(n):
        for j in range(i + 1, n):
            if abs(J[i, j]) > 1e-14:
                terms.append(("ZZ", [i, j], J[i, j]))
    # X terms
    for i in range(n):
        terms.append(("X", [i], h))
    return SparsePauliOp.from_sparse_list(terms, num_qubits=n)


def magnetization_op(n: int) -> SparsePauliOp:
    """Build  M_z = (1/N) sum_i Z_i."""
    terms = [("Z", [i], 1.0 / n) for i in range(n)]
    return SparsePauliOp.from_sparse_list(terms, num_qubits=n)


# ---------------------------------------------------------------------------
# Exact Trotterized time evolution  (for small N, used as reference)
# ---------------------------------------------------------------------------

def _op_to_matrix(op: SparsePauliOp) -> np.ndarray:
    """Convert SparsePauliOp to dense matrix."""
    return op.to_matrix(sparse=False)


def exact_trotter_magnetization(
    n: int,
    alpha: float,
    h: float,
    theta: float,
    n_steps: int,
    geometry: str = "chain",
    cutoff_frac: float = 1e-14,
) -> float:
    """Compute <M_z> after n_steps of first-order Trotter evolution.

    Each Trotter step applies:
        U(theta) = prod_{i<j} exp(-i * J_ij * Z_i Z_j)  *  prod_i exp(-i * theta/2 * X_i)

    followed by a final X-rotation layer (to match the IBM tutorial convention).

    Uses the same coupling cutoff as the circuit builders, so the exact
    reference matches the circuit being simulated.

    Parameters
    ----------
    n : int           Number of qubits.
    alpha : float     Power-law exponent.
    h : float         Field strength (unused in Trotter gate — theta controls kick).
    theta : float     Kick angle.
    n_steps : int     Number of Trotter steps.
    geometry : str    "chain" or "honeycomb".
    cutoff_frac : float  Coupling cutoff (fraction of J_max).  Must match
                         the value used in build_ideal_trotter_circuit.

    Returns
    -------
    float : expectation value of M_z.
    """
    dim = 2 ** n

    # Build ZZ unitary layer — using the same cutoff as the circuits
    pairs = coupling_pairs(n, alpha, cutoff_frac=cutoff_frac, geometry=geometry)
    H_zz = SparsePauliOp.from_sparse_list(
        [("ZZ", [i, j], J_ij) for i, j, J_ij in pairs],
        num_qubits=n,
    )
    U_zz = expm(-1j / 2 * _op_to_matrix(H_zz))  # Qiskit's Rzz(J) = exp(-i*J/2*ZZ)

    # Build X kick layer
    H_x = SparsePauliOp.from_sparse_list(
        [("X", [i], 1.0) for i in range(n)],
        num_qubits=n,
    )
    U_x = expm(-1j * (theta / 2.0) * _op_to_matrix(H_x))

    # Initial state |00...0>
    psi = np.zeros(dim, dtype=complex)
    psi[0] = 1.0

    # Apply Trotter steps
    for _ in range(n_steps):
        psi = U_x @ psi
        psi = U_zz @ psi

    # Final X rotation
    psi = U_x @ psi

    # Measure <M_z>
    Mz = _op_to_matrix(magnetization_op(n))
    return np.real(psi.conj() @ Mz @ psi)


def exact_magnetization_sweep(
    n: int,
    alpha: float,
    h: float,
    thetas: np.ndarray,
    max_steps: int,
    geometry: str = "chain",
    cutoff_frac: float = 0.10,
) -> np.ndarray:
    """Compute exact <M_z> for a grid of (theta, trotter_step).

    Returns
    -------
    result : np.ndarray of shape (len(thetas), max_steps + 1)
        result[t, s] = <M_z> at theta=thetas[t] after s Trotter steps.
    """
    result = np.zeros((len(thetas), max_steps + 1))
    for t_idx, theta in enumerate(thetas):
        for s in range(max_steps + 1):
            result[t_idx, s] = exact_trotter_magnetization(
                n, alpha, h, theta, s, geometry=geometry,
                cutoff_frac=cutoff_frac,
            )
    return result
