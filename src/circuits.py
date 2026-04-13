"""
Connectivity-tax benchmark — circuit construction.

Three variants of Trotterized time-evolution circuits:
1. Ideal (all-to-all connectivity)
2. Heavy-hex unitary (transpiled with SWAP routing)
3. Heavy-hex dynamic (ancilla + mid-circuit measurement + feedforward)
"""

import numpy as np
from qiskit.circuit import QuantumCircuit, QuantumRegister, ClassicalRegister, Parameter
from qiskit.circuit.classical import expr
from qiskit.quantum_info import SparsePauliOp
from qiskit.transpiler import CouplingMap
from qiskit.transpiler.preset_passmanagers import generate_preset_pass_manager

from .model import coupling_matrix, coupling_pairs, distance_matrix


# ---------------------------------------------------------------------------
# Optimal gate scheduling via edge coloring
# ---------------------------------------------------------------------------

def _optimal_rzz_schedule(pairs):
    """Find minimum-depth parallel schedule for Rzz gates.

    Groups coupling pairs into layers where no two pairs in the same
    layer share a qubit.  Uses exact edge coloring (backtracking) to
    find the minimum number of layers, which equals the chromatic index
    of the interaction graph.

    All Rzz gates commute ([Z_i Z_j, Z_k Z_l] = 0), so reordering
    is always safe.

    Parameters
    ----------
    pairs : list of (i, j, J_ij)

    Returns
    -------
    list of layers, each layer is a list of (i, j, J_ij).
    """
    if not pairs:
        return []

    import rustworkx as rx

    # Build graph: nodes = qubits, edges = pairs
    all_qubits = sorted({q for i, j, _ in pairs for q in (i, j)})
    q_to_idx = {q: idx for idx, q in enumerate(all_qubits)}

    g = rx.PyGraph()
    for q in all_qubits:
        g.add_node(q)

    edge_info = {}
    for i, j, J in pairs:
        eidx = g.add_edge(q_to_idx[i], q_to_idx[j], None)
        edge_info[eidx] = (i, j, J)

    edge_list = list(g.edge_list())
    n_edges = len(edge_list)
    max_deg = max(g.degree(n) for n in g.node_indices())

    # Backtracking edge coloring: try max_degree first, then max_degree + 1
    def try_color(n_colors):
        coloring = [-1] * n_edges

        def conflicts(eidx, color):
            src, tgt = edge_list[eidx]
            for other in range(n_edges):
                if other == eidx or coloring[other] != color:
                    continue
                os, ot = edge_list[other]
                if src in (os, ot) or tgt in (os, ot):
                    return True
            return False

        def solve(idx):
            if idx == n_edges:
                return True
            for c in range(n_colors):
                if not conflicts(idx, c):
                    coloring[idx] = c
                    if solve(idx + 1):
                        return True
                    coloring[idx] = -1
            return False

        return coloring if solve(0) else None

    for n_colors in range(max_deg, max_deg + 2):
        result = try_color(n_colors)
        if result is not None:
            layers = [[] for _ in range(n_colors)]
            for eidx, color in enumerate(result):
                layers[color].append(edge_info[eidx])
            return [l for l in layers if l]

    # Fallback to Misra-Gries (should never reach here for small graphs)
    edge_colors = rx.graph_misra_gries_edge_color(g)
    n_c = max(edge_colors.values()) + 1
    layers = [[] for _ in range(n_c)]
    for eidx, color in edge_colors.items():
        layers[color].append(edge_info[eidx])
    return [l for l in layers if l]


# ---------------------------------------------------------------------------
# Variant 1: Ideal circuit (all-to-all connectivity)
# ---------------------------------------------------------------------------

def build_ideal_trotter_circuit(
    n: int,
    alpha: float,
    theta: Parameter | float,
    n_steps: int = 1,
    measure: bool = True,
    final_rot: bool = True,
    cutoff_frac: float = 0.10,
    geometry: str = "chain",
) -> QuantumCircuit:
    """Build Trotter circuit assuming all-to-all connectivity.

    Each step:  Rx(theta) on all qubits, then Rzz(J_ij) on all coupled pairs.
    Rzz gates are scheduled in optimal parallel layers via edge coloring.
    Final Rx(theta) layer appended if final_rot=True.

    Parameters
    ----------
    n : int             Number of qubits.
    alpha : float       Power-law exponent.
    theta :             Kick angle (Parameter or float).
    n_steps : int       Number of Trotter steps.
    measure : bool      Append measurements.
    final_rot : bool    Append final Rx layer.
    cutoff_frac : float Coupling cutoff (fraction of J_max).
    geometry : str      "chain" or "honeycomb".

    Returns
    -------
    QuantumCircuit
    """
    pairs = coupling_pairs(n, alpha, cutoff_frac=cutoff_frac, geometry=geometry)
    layers = _optimal_rzz_schedule(pairs)
    qc = QuantumCircuit(n)

    for _ in range(n_steps):
        # X kick
        for i in range(n):
            qc.rx(theta, i)
        qc.barrier()
        # ZZ interactions in optimal parallel layers
        for layer in layers:
            for (i, j, J_ij) in layer:
                qc.rzz(J_ij, i, j)
            qc.barrier()

    if final_rot:
        for i in range(n):
            qc.rx(theta, i)

    if measure:
        qc.measure_all()
    return qc


# ---------------------------------------------------------------------------
# Connectivity helpers
# ---------------------------------------------------------------------------

def ideal_coupling_map(n: int) -> CouplingMap:
    """All-to-all coupling map for n qubits."""
    edges = []
    for i in range(n):
        for j in range(i + 1, n):
            edges.append((i, j))
            edges.append((j, i))
    return CouplingMap(edges)


def heavyhex_coupling_map_small() -> CouplingMap:
    """Return a small heavy-hex coupling map (IBM Eagle / Heron style).

    This is a manually specified 27-qubit heavy-hex fragment that matches
    the topology of ibm_brisbane / ibm_osaka style backends.  For our
    4-8 qubit experiments we only use a subgraph, but the transpiler
    needs the full map.

    Qubit layout (heavy-hex pattern):
        Degree-3 nodes: 0,1,2,3,4,5,6,7,8,9,10,11,12 (hex lattice sites)
        Degree-2 nodes: 13,14,15,16,17,18,19,20,21,22,23,24,25,26 (bridge qubits)

    For simplicity, we use CouplingMap.from_heavy_hex() when available,
    otherwise build manually.
    """
    try:
        return CouplingMap.from_heavy_hex(3, bidirectional=True)
    except AttributeError:
        pass

    # Manual heavy-hex for 27 qubits (3x3 heavy hex)
    # This is the connectivity of ibm_hanoi / ibm_cairo style
    edges = [
        (0, 1), (1, 2), (2, 3), (3, 4),
        (0, 14), (4, 18),
        (14, 9), (15, 9), (15, 5), (16, 5), (16, 10),
        (17, 10), (17, 6), (18, 6),
        (9, 20), (10, 23),
        (20, 11), (21, 11), (21, 7), (22, 7), (22, 12),
        (23, 12), (23, 8), (24, 8),
        (1, 13), (13, 9), (2, 15), (3, 17),
        (5, 19), (19, 11), (6, 21),
        (7, 25), (25, 12), (8, 26),
    ]
    bidirectional = []
    for (a, b) in edges:
        bidirectional.append((a, b))
        bidirectional.append((b, a))
    return CouplingMap(bidirectional)


def linear_coupling_map(n: int) -> CouplingMap:
    """Simple linear chain coupling map."""
    edges = []
    for i in range(n - 1):
        edges.append((i, i + 1))
        edges.append((i + 1, i))
    return CouplingMap(edges)


# ---------------------------------------------------------------------------
# Variant 2: Transpile to heavy-hex with SWAP routing
# ---------------------------------------------------------------------------

def transpile_to_heavyhex(
    circuit: QuantumCircuit,
    coupling_map: CouplingMap | None = None,
    optimization_level: int = 1,
    seed: int = 42,
    n_seeds: int = 1,
) -> QuantumCircuit:
    """Transpile a circuit to a heavy-hex coupling map.

    Parameters
    ----------
    circuit : QuantumCircuit
        The logical circuit (e.g. from build_ideal_trotter_circuit).
    coupling_map : CouplingMap or None
        Target topology.  If None, uses a default heavy-hex map.
    optimization_level : int
        Qiskit transpiler optimization level (0-3).
    seed : int
        Transpiler seed for reproducibility (used when n_seeds=1).
    n_seeds : int
        If > 1, try this many random seeds and return the circuit with
        the lowest two-qubit depth.  Uses optimization_level=3 for
        multi-seed runs.

    Returns
    -------
    QuantumCircuit : the transpiled (physical) circuit.
    """
    if coupling_map is None:
        coupling_map = heavyhex_coupling_map_small()

    if n_seeds <= 1:
        pm = generate_preset_pass_manager(
            optimization_level=optimization_level,
            coupling_map=coupling_map,
            seed_transpiler=seed,
        )
        return pm.run(circuit)

    # Multi-seed: try many seeds, keep the best by 2Q depth
    two_q_set = {"cx", "cz", "ecr", "rzz", "rxx", "ryy", "cp", "swap", "iswap"}
    opt = max(optimization_level, 3)  # use highest optimization for multi-seed
    best_circ = None
    best_depth = float("inf")

    for s in range(n_seeds):
        pm = generate_preset_pass_manager(
            optimization_level=opt,
            coupling_map=coupling_map,
            seed_transpiler=s,
        )
        candidate = pm.run(circuit)

        # Count 2Q depth
        from qiskit.converters import circuit_to_dag
        dag = circuit_to_dag(candidate)
        depth = sum(
            1 for layer in dag.layers()
            if any(
                len(node.qargs) >= 2 and node.op.name.lower() in two_q_set
                for node in layer["graph"].op_nodes()
            )
        )
        if depth < best_depth:
            best_depth = depth
            best_circ = candidate

    return best_circ


def transpile_to_ideal(
    circuit: QuantumCircuit,
    n: int,
    optimization_level: int = 1,
    seed: int = 42,
) -> QuantumCircuit:
    """Transpile a circuit with all-to-all connectivity (baseline)."""
    pm = generate_preset_pass_manager(
        optimization_level=optimization_level,
        coupling_map=ideal_coupling_map(n),
        seed_transpiler=seed,
    )
    return pm.run(circuit)


# ---------------------------------------------------------------------------
# Variant 3: Dynamic circuit — ancilla-based ZZ via measure + feedforward
# ---------------------------------------------------------------------------

def _color_pairs(pairs):
    """Deprecated: use _optimal_rzz_schedule instead."""
    return _optimal_rzz_schedule(pairs)


def _build_dynamic_rzz_block(
    qc: QuantumCircuit,
    qi: int,
    qj: int,
    ancilla: int,
    creg_bit,
    angle: float,
):
    """Append a dynamic-circuit Rzz(angle) between qi and qj using an ancilla.

    Protocol:
    1. Reset ancilla to |0>
    2. CX qi -> ancilla, CX qj -> ancilla  (parity computation)
    3. Rz(angle) on ancilla
    4. H + Measure ancilla  (X-basis measurement)
    5. If result == 1: apply Z to qi and Z to qj  (correction)

    This effectively implements  Rzz(angle) = exp(-i * angle/2 * Z_i Z_j).
    """
    qc.reset(ancilla)
    qc.cx(qi, ancilla)
    qc.cx(qj, ancilla)
    qc.rz(angle, ancilla)  # Rz(J) gives effective exp(-i*J/2*ZZ), matching Rzz(J)
    qc.h(ancilla)
    qc.measure(ancilla, creg_bit)
    with qc.if_test((creg_bit, 1)):
        qc.z(qi)
        qc.z(qj)


def build_dynamic_trotter_circuit(
    n: int,
    alpha: float,
    theta: Parameter | float,
    n_steps: int = 1,
    measure: bool = True,
    final_rot: bool = True,
    cutoff_frac: float = 0.10,
    geometry: str = "chain",
) -> QuantumCircuit:
    """Build hybrid Trotter circuit on *logical* qubits (for noiseless
    simulation and reference values).

    NN pairs use the ancilla-assisted protocol; non-NN pairs use bare Rzz.
    Each NN pair gets a dedicated ancilla.

    For transpilation to a specific coupling map, use
    ``build_dynamic_trotter_circuit_physical`` instead — it places NN
    ancilla blocks directly on physical bridge qubits so the transpiler
    cannot insert SWAPs for them.

    Parameters
    ----------
    n : int             Number of data qubits.
    alpha : float       Power-law exponent.
    theta :             Kick angle.
    n_steps : int       Trotter steps.
    measure : bool      Final measurement of data qubits.
    final_rot : bool    Final Rx layer.
    cutoff_frac : float Coupling cutoff (fraction of J_max).
    geometry : str      "chain" or "honeycomb".

    Returns
    -------
    QuantumCircuit with n data + len(nn_pairs) ancilla qubits.
    """
    pairs = coupling_pairs(n, alpha, cutoff_frac=cutoff_frac, geometry=geometry)

    # Split pairs into nearest-neighbour (distance == 1) and long-range
    D = distance_matrix(n, geometry)
    nn_pairs = [(i, j, J) for i, j, J in pairs if D[i, j] == 1.0]
    lr_pairs = [(i, j, J) for i, j, J in pairs if D[i, j] > 1.0]

    nn_layers = _optimal_rzz_schedule(nn_pairs) if nn_pairs else []
    lr_layers = _optimal_rzz_schedule(lr_pairs) if lr_pairs else []

    # One dedicated ancilla per NN pair
    n_anc = len(nn_pairs)
    pair_to_anc = {(i, j): k for k, (i, j, _) in enumerate(nn_pairs)}

    data_reg = QuantumRegister(n, name="d")
    qc_regs = [data_reg]
    if n_anc > 0:
        anc_reg = QuantumRegister(n_anc, name="anc")
        mid_creg = ClassicalRegister(n_anc, name="mid")
        qc_regs += [anc_reg, mid_creg]
    else:
        anc_reg = None
        mid_creg = None

    qc = QuantumCircuit(*qc_regs)

    for _ in range(n_steps):
        for i in range(n):
            qc.rx(theta, data_reg[i])
        qc.barrier()

        for layer in nn_layers:
            for (i, j, J_ij) in layer:
                k = pair_to_anc[(i, j)]
                _build_dynamic_rzz_block(
                    qc, data_reg[i], data_reg[j],
                    anc_reg[k], mid_creg[k], J_ij,
                )
            qc.barrier()

        for layer in lr_layers:
            for (i, j, J_ij) in layer:
                qc.rzz(J_ij, data_reg[i], data_reg[j])
            qc.barrier()

    if final_rot:
        for i in range(n):
            qc.rx(theta, data_reg[i])

    if measure:
        final_creg = ClassicalRegister(n, name="meas")
        qc.add_register(final_creg)
        for i in range(n):
            qc.measure(data_reg[i], final_creg[i])

    return qc


def build_dynamic_trotter_circuit_physical(
    n: int,
    alpha: float,
    theta: Parameter | float,
    coupling_map: CouplingMap,
    n_steps: int = 1,
    measure: bool = True,
    final_rot: bool = True,
    cutoff_frac: float = 0.10,
    geometry: str = "chain",
) -> QuantumCircuit:
    """Build hybrid Trotter circuit directly on physical qubits of a
    heavy-hex coupling map.

    The circuit register spans the **full** coupling map.  Data qubits
    sit on degree-3 nodes, and each NN ancilla sits on the degree-2
    bridge qubit between its two data qubits.  Because every CX in the
    ancilla-assisted ZZ block targets qubits that are already adjacent
    on the hardware graph, the transpiler's routing pass will never
    insert SWAPs for NN interactions.

    Non-NN pairs are emitted as bare Rzz gates between (non-adjacent)
    data-node qubits; the transpiler routes those with SWAPs, exactly
    as the unitary variant does.

    Use ``transpile_dynamic_physical`` to transpile the returned circuit
    (basis-gate decomposition + routing of LR gates only).

    Parameters
    ----------
    n : int               Number of data qubits.
    alpha : float         Power-law exponent.
    theta :               Kick angle.
    coupling_map :        Target heavy-hex CouplingMap.
    n_steps : int         Trotter steps.
    measure : bool        Final measurement of data qubits.
    final_rot : bool      Final Rx layer.
    cutoff_frac : float   Coupling cutoff (fraction of J_max).
    geometry : str        "chain" or "honeycomb".

    Returns
    -------
    QuantumCircuit on ``coupling_map.size()`` physical qubits.
    ``circuit.metadata["data_physical_qubits"]`` lists the physical
    qubit indices that carry the data (needed for observable measurement).
    """
    # ---- physical layout on the coupling map ----
    if geometry == "honeycomb" and n == 6:
        data_nodes, bridge_list = find_heavyhex_hexring(coupling_map)
    else:
        data_nodes, bridge_list = find_heavyhex_chain(coupling_map, n)

    # bridge_lookup: (logical_i, logical_j) -> physical bridge qubit
    bridge_lookup = {}
    for k in range(len(bridge_list)):
        if geometry == "honeycomb" and n == 6:
            a, b = k, (k + 1) % n
        else:
            a, b = k, k + 1
        bridge_lookup[(a, b)] = bridge_list[k]
        bridge_lookup[(b, a)] = bridge_list[k]

    # ---- split coupling pairs ----
    pairs = coupling_pairs(n, alpha, cutoff_frac=cutoff_frac, geometry=geometry)
    D = distance_matrix(n, geometry)
    nn_pairs = [(i, j, J) for i, j, J in pairs if D[i, j] == 1.0]
    lr_pairs = [(i, j, J) for i, j, J in pairs if D[i, j] > 1.0]

    nn_layers = _optimal_rzz_schedule(nn_pairs) if nn_pairs else []
    lr_layers = _optimal_rzz_schedule(lr_pairs) if lr_pairs else []

    # ---- build circuit on the full physical topology ----
    n_phys = coupling_map.size()
    qr = QuantumRegister(n_phys, name="q")

    n_nn = len(nn_pairs)
    pair_to_mid = {(i, j): k for k, (i, j, _) in enumerate(nn_pairs)}

    regs = [qr]
    mid_creg = None
    if n_nn > 0:
        mid_creg = ClassicalRegister(n_nn, name="mid")
        regs.append(mid_creg)

    qc = QuantumCircuit(*regs)

    for _ in range(n_steps):
        # X kick on data qubits (physical indices)
        for i in range(n):
            qc.rx(theta, qr[data_nodes[i]])
        qc.barrier()

        # NN interactions — ancilla blocks on physical bridge qubits
        for layer in nn_layers:
            for (i, j, J_ij) in layer:
                p_i = data_nodes[i]
                p_j = data_nodes[j]
                p_anc = bridge_lookup[(i, j)]
                mb = pair_to_mid[(i, j)]
                # ancilla-assisted Rzz directly on adjacent hardware qubits
                qc.reset(qr[p_anc])
                qc.cx(qr[p_i], qr[p_anc])
                qc.cx(qr[p_j], qr[p_anc])
                qc.rz(J_ij, qr[p_anc])
                qc.h(qr[p_anc])
                qc.measure(qr[p_anc], mid_creg[mb])
                with qc.if_test((mid_creg[mb], 1)):
                    qc.z(qr[p_i])
                    qc.z(qr[p_j])
            qc.barrier()

        # Long-range interactions — bare Rzz, transpiler routes these
        for layer in lr_layers:
            for (i, j, J_ij) in layer:
                qc.rzz(J_ij, qr[data_nodes[i]], qr[data_nodes[j]])
            qc.barrier()

    if final_rot:
        for i in range(n):
            qc.rx(theta, qr[data_nodes[i]])

    if measure:
        meas_creg = ClassicalRegister(n, name="meas")
        qc.add_register(meas_creg)
        for i in range(n):
            qc.measure(qr[data_nodes[i]], meas_creg[i])

    qc.metadata = {
        "data_physical_qubits": [data_nodes[i] for i in range(n)],
        "bridge_physical_qubits": [bridge_lookup[(i, j)]
                                   for i, j, _ in nn_pairs],
        "nn_pairs_order": [(i, j) for i, j, _ in nn_pairs],
        "geometry": geometry,
    }

    return qc


def transpile_dynamic_physical(
    circuit: QuantumCircuit,
    coupling_map: CouplingMap | None = None,
    optimization_level: int = 1,
    seed: int = 42,
    backend=None,
) -> QuantumCircuit:
    """Transpile a physical-qubit dynamic circuit.

    Because the circuit was built by ``build_dynamic_trotter_circuit_physical``
    directly on the coupling-map topology, the identity layout is correct
    and the routing pass only needs to handle the long-range Rzz gates
    (if any).  NN ancilla blocks are already between adjacent physical
    qubits and pass through routing untouched.

    Parameters
    ----------
    circuit : QuantumCircuit
        Output of ``build_dynamic_trotter_circuit_physical``.
    coupling_map : CouplingMap or None
        Same coupling map used to build the circuit.  Ignored when
        *backend* is provided.
    optimization_level : int
        Qiskit transpiler optimization level (0-3).
    seed : int
        Transpiler seed for reproducibility.
    backend : optional
        A Qiskit backend (real or fake).  When provided, the pass
        manager uses the backend's full target — basis gates, timing,
        calibrations — so the output circuit is ISA-compliant and can
        be submitted directly to the Sampler.  If None, falls back to
        *coupling_map* (sufficient for local simulation / metric
        extraction, but the output may contain non-native gates).

    Returns
    -------
    QuantumCircuit : transpiled circuit.
    """
    trivial_layout = list(range(circuit.num_qubits))

    pm_kwargs = dict(
        optimization_level=optimization_level,
        initial_layout=trivial_layout,
        seed_transpiler=seed,
    )
    if backend is not None:
        pm_kwargs["backend"] = backend
    elif coupling_map is not None:
        pm_kwargs["coupling_map"] = coupling_map
    else:
        raise ValueError("Provide either coupling_map or backend.")

    pm = generate_preset_pass_manager(**pm_kwargs)
    result = pm.run(circuit)

    # Propagate metadata so simulate_mz can find data qubits
    if circuit.metadata:
        if result.metadata is None:
            result.metadata = {}
        result.metadata.update(circuit.metadata)

    return result


def magnetization_observable(n: int) -> SparsePauliOp:
    """M_z = (1/N) sum_i Z_i  for the data qubits."""
    terms = [("Z", [i], 1.0 / n) for i in range(n)]
    return SparsePauliOp.from_sparse_list(terms, num_qubits=n)


def magnetization_observable_dynamic(n_data: int, n_anc: int) -> SparsePauliOp:
    """M_z for dynamic circuit: n_data data qubits + n_anc ancillas."""
    total = n_data + n_anc
    terms = [("Z", [i], 1.0 / n_data) for i in range(n_data)]
    return SparsePauliOp.from_sparse_list(terms, num_qubits=total)


# ---------------------------------------------------------------------------
# Heavy-hex embedding helpers for manual initial layout
# ---------------------------------------------------------------------------

def _get_heavyhex_bridges(coupling_map: CouplingMap):
    """Extract degree-3 adjacency and bridge map from a heavy-hex coupling map.

    Returns
    -------
    deg3 : list of int       Degree-3 (data) nodes.
    bridges : dict           (sorted pair of deg3 nodes) -> deg2 bridge qubit.
    adj : dict               deg3 node -> list of deg3 neighbours.
    """
    import rustworkx as rx
    g = coupling_map.graph.to_undirected(multigraph=False)
    degrees = {n: g.degree(n) for n in g.node_indices()}
    deg3 = sorted(n for n, d in degrees.items() if d == 3)
    deg2 = sorted(n for n, d in degrees.items() if d == 2)

    bridges = {}
    adj = {n: [] for n in deg3}
    for d2 in deg2:
        neighbors = [n for n in g.neighbors(d2) if degrees.get(n, 0) == 3]
        if len(neighbors) == 2:
            a, b = sorted(neighbors)
            bridges[(a, b)] = d2
            adj[a].append(b)
            adj[b].append(a)
    return deg3, bridges, adj


def find_heavyhex_hexring(coupling_map: CouplingMap):
    """Find a hexagonal ring of 6 degree-3 nodes on a heavy-hex map.

    Returns
    -------
    ring : list of 6 int
        Physical qubit indices forming a hexagonal ring.
    ring_bridges : list of 6 int
        Physical qubit indices for the bridge (ancilla) between
        consecutive ring nodes.  ring_bridges[i] sits between
        ring[i] and ring[(i+1)%6].
    """
    deg3, bridges, adj = _get_heavyhex_bridges(coupling_map)

    def find_ring(start, path=None, visited=None):
        if path is None:
            path = [start]
            visited = {start}
        if len(path) == 6:
            if start in adj[path[-1]]:
                return path
            return None
        for nb in adj[path[-1]]:
            if nb not in visited:
                result = find_ring(start, path + [nb], visited | {nb})
                if result:
                    return result
        return None

    for s in deg3:
        ring = find_ring(s)
        if ring:
            ring_bridges = []
            for i in range(6):
                key = tuple(sorted([ring[i], ring[(i + 1) % 6]]))
                ring_bridges.append(bridges[key])
            return ring, ring_bridges

    raise ValueError("Cannot find a hexagonal ring on this heavy-hex map")


def find_heavyhex_chain(coupling_map: CouplingMap, n_data: int):
    """Find an optimal chain of n_data degree-3 nodes on a heavy-hex map.

    Returns
    -------
    data_chain : list of int
        Physical qubit indices for data qubits (degree-3 nodes).
    nn_bridges : list of int
        Physical qubit indices for NN ancillas (degree-2 bridge between
        consecutive data qubits).  Length = n_data - 1.
    """
    g = coupling_map.graph.to_undirected(multigraph=False)
    degrees = {n: g.degree(n) for n in g.node_indices()}
    deg3 = sorted(n for n, d in degrees.items() if d == 3)
    deg2 = sorted(n for n, d in degrees.items() if d == 2)

    # Build bridge map: (deg3_a, deg3_b) -> deg2_bridge
    bridges = {}
    for d2 in deg2:
        neighbors = [n for n in g.neighbors(d2) if degrees.get(n, 0) == 3]
        if len(neighbors) == 2:
            a, b = sorted(neighbors)
            bridges[(a, b)] = d2

    # Adjacency among deg-3 nodes
    adj = {n: [] for n in deg3}
    for (a, b) in bridges:
        adj[a].append(b)
        adj[b].append(a)

    # DFS for a chain of the required length
    def find_chain(start, length, visited=None):
        if visited is None:
            visited = [start]
        if len(visited) == length:
            return visited
        for nb in adj[start]:
            if nb not in visited:
                result = find_chain(nb, length, visited + [nb])
                if result:
                    return result
        return None

    for s in deg3:
        chain = find_chain(s, n_data)
        if chain:
            bridge_list = []
            for i in range(len(chain) - 1):
                key = tuple(sorted([chain[i], chain[i + 1]]))
                bridge_list.append(bridges[key])
            return chain, bridge_list

    raise ValueError(f"Cannot find a chain of {n_data} degree-3 nodes on this map")


def transpile_dynamic_with_layout(
    circuit: QuantumCircuit,
    coupling_map: CouplingMap,
    n_data: int,
    geometry: str = "chain",
    optimization_level: int = 1,
    seed: int = 42,
) -> QuantumCircuit:
    """Transpile a dynamic circuit with manual initial layout on heavy-hex.

    Uses ``circuit.metadata["nn_pairs_order"]`` (set by
    ``build_dynamic_trotter_circuit``) to pin each ancilla to the hardware
    bridge qubit that sits between its two data qubits.  This guarantees
    that all CX gates in the ancilla-assisted ZZ blocks are already
    nearest-neighbour on the hardware, so the transpiler inserts **zero**
    SWAPs for the NN interactions.

    For geometry="chain": places data qubits on a chain of degree-3 nodes.
    For geometry="honeycomb" (n_data=6): places data qubits on a hexagonal
    ring of degree-3 nodes.
    """
    n_anc = circuit.num_qubits - n_data

    if geometry == "honeycomb" and n_data == 6:
        ring, ring_bridges = find_heavyhex_hexring(coupling_map)
        data_nodes = ring
        bridge_nodes = ring_bridges
    else:
        data_nodes, bridge_nodes = find_heavyhex_chain(coupling_map, n_data)

    # Build lookup: (logical_i, logical_j) -> physical bridge qubit
    # For chain: bridge_nodes[k] sits between data_nodes[k] and data_nodes[k+1]
    # For hex ring: bridge_nodes[k] sits between data_nodes[k] and data_nodes[(k+1)%6]
    bridge_lookup = {}
    n_bridges = len(bridge_nodes)
    for k in range(n_bridges):
        if geometry == "honeycomb" and n_data == 6:
            a, b = k, (k + 1) % n_data
        else:
            a, b = k, k + 1
        # Store both orderings for easy lookup
        bridge_lookup[(a, b)] = bridge_nodes[k]
        bridge_lookup[(b, a)] = bridge_nodes[k]

    # Build initial layout
    initial_layout = list(range(circuit.num_qubits))

    # Map data qubits: logical qubit i -> data_nodes[i]
    for i in range(n_data):
        initial_layout[i] = data_nodes[i]

    # Map ancillas using the nn_pairs_order metadata
    nn_pairs_order = []
    if circuit.metadata and "nn_pairs_order" in circuit.metadata:
        nn_pairs_order = circuit.metadata["nn_pairs_order"]

    used = set(data_nodes)
    anc_positions = []

    for k, (li, lj) in enumerate(nn_pairs_order):
        # Pin ancilla k to the bridge between logical qubits li and lj
        bridge = bridge_lookup.get((li, lj))
        if bridge is not None and bridge not in used:
            anc_positions.append(bridge)
            used.add(bridge)
        else:
            # Fallback: shouldn't happen for well-formed NN pairs
            anc_positions.append(None)

    # Fill any remaining ancillas (from fallback or extra ancillas)
    # with unused physical qubits
    all_physical = set(range(coupling_map.size()))
    remaining = sorted(all_physical - used)
    for k in range(len(anc_positions)):
        if anc_positions[k] is None:
            if remaining:
                q = remaining.pop(0)
                anc_positions[k] = q
                used.add(q)
    while len(anc_positions) < n_anc:
        if remaining:
            q = remaining.pop(0)
            anc_positions.append(q)
            used.add(q)
        else:
            break

    for i in range(n_anc):
        initial_layout[n_data + i] = anc_positions[i]

    pm = generate_preset_pass_manager(
        optimization_level=optimization_level,
        coupling_map=coupling_map,
        initial_layout=initial_layout,
        seed_transpiler=seed,
    )
    return pm.run(circuit)
