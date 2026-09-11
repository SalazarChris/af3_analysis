"""
V3 Automatic Site/Region Detection.

Data-driven definitions for the local-geometry (F10) and domain-motion (F11)
figures, so that sites and regions do not have to be supplied manually.

Methodology (descriptive/structural only; no functional annotation):

Sites (displacement-concentration detection):
    1. For every seed, compute per-residue CA displacement between each
       target condition and the reference condition (matched seed/sample
       pairing, identical to figure F04).
    2. Score each residue by its mean displacement across seeds (per-seed
       means first, so every seed contributes equally regardless of the
       number of samples it contains).
    3. Flag residues whose score is a robust outlier (z >= threshold using
       median/MAD; degenerate all-tied distributions flag every score that
       differs from the bulk level).
    4. Group flagged residues within a chain when they are within
       `grouping_gap` sequence positions of each other. Each group becomes
       one site: center = highest-scoring residue, radius = the configured
       contact distance, label = site_N ordered by descending peak score.

Regions (contact-graph component detection):
    1. Build the residue contact graph of the reference structure
       (CA-CA distance < contact threshold), excluding sequence-adjacent
       residue pairs.
    2. Connected components with at least `min_residues` residues become
       regions, reported as seq_list plus start/end for reporting.

All outputs are *predicted structural sites/regions* derived from the data.
Callers must surface this in warnings/manifests; detection failure returns
an empty list so figures skip rather than inventing definitions.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from ..model import StructureData
from .displacement import calculate_per_residue_displacement


def detect_displacement_sites(
    dataset: Any,
    ref_condition: Optional[str],
    *,
    z_threshold: float = 2.0,
    radius: float = 8.0,
    alignment_atom: str = "CA",
    grouping_gap: int = 2,
    max_sites: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """
    Detect sites of concentrated between-condition displacement.

    Parameters
    ----------
    dataset : Dataset
        V3 normalized dataset (condition -> seed -> sample -> structure).
    ref_condition : str or None
        Reference condition id. All target conditions are compared to it.
    z_threshold : float
        Robust z-score threshold for a residue to be a displacement peak
        (deliberately conservative-positive: this is exploratory flagging
        for further inspection, not a significance test).
    radius : float
        Radius (Angstrom) attached to each detected site (uses the
        pipeline's configured contact distance by default).
    alignment_atom : str
        Atom used for displacement (passed through unchanged).
    grouping_gap : int
        Flagged residues within this many sequence positions of each other
        (same chain) are grouped into a single site.
    max_sites : int, optional
        Keep only the top-N sites by peak displacement score.

    Returns
    -------
    list of dict compatible with calculate_local_geometry region definitions:
        {"label", "chain", "residue", "radius", "mean_displacement", "z"}
    Empty if no reference condition, no target conditions, no matched
    predictions, or no significant displacement concentration.
    """
    if ref_condition is None:
        return []

    predictions = getattr(dataset, "predictions", {}) or {}
    target_conditions = sorted(
        cid for cid in predictions.keys() if cid != ref_condition
    )
    if not target_conditions:
        return []

    # Per-seed displacement fields: seed_key -> residue_key -> [displacements]
    seed_fields: Dict[Any, Dict[Any, List[float]]] = {}

    for condition_id in target_conditions:
        for seed in sorted(predictions.get(condition_id, {}).keys()):
            seed_key = (condition_id, seed)
            for sample in sorted(predictions[condition_id][seed].keys()):
                target_struct = predictions[condition_id][seed][sample]
                ref_struct = (
                    predictions.get(ref_condition, {})
                    .get(seed, {})
                    .get(sample)
                )
                if ref_struct is None or target_struct is None:
                    continue
                try:
                    disp = calculate_per_residue_displacement(
                        ref_struct,
                        target_struct,
                        alignment_atom=alignment_atom,
                    )
                except Exception:
                    continue
                field = seed_fields.setdefault(seed_key, {})
                for row in disp["displacements"]:
                    key = (row["chain_id"], row["residue_index"])
                    field.setdefault(key, []).append(float(row["displacement"]))

    if not seed_fields:
        return []

    # Per-residue score: mean across seeds of per-seed mean displacement.
    residue_scores: Dict[Any, float] = {}
    all_keys = set()
    for field in seed_fields.values():
        all_keys.update(field.keys())
    for key in all_keys:
        seed_means = [
            float(np.mean(field[key]))
            for field in seed_fields.values()
            if key in field
        ]
        if seed_means:
            residue_scores[key] = float(np.mean(seed_means))

    if not residue_scores:
        return []

    keys = sorted(residue_scores.keys())
    scores = np.array([residue_scores[k] for k in keys], dtype=np.float64)

    flagged = _robust_outlier_indices(scores, z_threshold)
    if not flagged:
        return []

    z_scores = _robust_z_scores(scores)

    # Group flagged residues per chain by sequence proximity.
    sites: List[Dict[str, Any]] = []
    for chain_id in sorted({keys[i][0] for i in flagged}):
        chain_flagged = sorted(
            (keys[i][1], i) for i in flagged if keys[i][0] == chain_id
        )
        group: List[tuple] = []
        groups: List[List[tuple]] = []
        previous_seq = None
        for seq_id, idx in chain_flagged:
            if (
                previous_seq is not None
                and seq_id - previous_seq > grouping_gap
            ):
                groups.append(group)
                group = []
            group.append((seq_id, idx))
            previous_seq = seq_id
        if group:
            groups.append(group)

        for g in groups:
            center_seq, center_idx = max(
                g, key=lambda t: (scores[t[1]], -t[0])
            )
            z_value = z_scores[center_idx]
            sites.append({
                "label": f"site_{len(sites) + 1}",
                "chain": chain_id,
                "residue": int(center_seq),
                "radius": float(radius),
                "mean_displacement": float(scores[center_idx]),
                "z": (float(z_value) if np.isfinite(z_value) else None),
            })

    # Relabel by descending peak displacement for stable ordering.
    sites.sort(key=lambda s: -s["mean_displacement"])
    for k, site in enumerate(sites, start=1):
        site["label"] = f"site_{k}"

    if max_sites is not None:
        sites = sites[:max_sites]

    return sites


def detect_contact_regions(
    reference_structure: StructureData,
    *,
    contact_threshold: float = 8.0,
    min_residues: int = 3,
) -> List[Dict[str, Any]]:
    """
    Detect regions as connected components of the reference contact graph.

    Parameters
    ----------
    reference_structure : StructureData
    contact_threshold : float
        CA-CA distance threshold defining a contact (Angstrom).
    min_residues : int
        Minimum number of residues for a component to be reported.

    Returns
    -------
    list of dict compatible with calculate_local_geometry region definitions:
        {"label", "chain", "start", "end", "seq_list", "n_residues"}
    Empty if the structure has no protein chains with CA coordinates or no
    sufficiently large contact-linked components.
    """
    regions: List[Dict[str, Any]] = []

    try:
        protein_chain_ids = reference_structure.get_protein_chains()
    except Exception:
        protein_chain_ids = []

    for chain_id in sorted(protein_chain_ids):
        chain = reference_structure.get_chain(chain_id)
        if chain is None:
            continue

        # Residues with CA coordinates, sorted by auth_seq_id.
        residues = sorted(
            (
                r.auth_seq_id
                for r in chain.residues
                if r.auth_seq_id is not None and r.ca_coords is not None
            )
        )
        if len(residues) < min_residues:
            continue

        coords = {
            r.auth_seq_id: np.array(r.ca_coords, dtype=np.float64)
            for r in chain.residues
            if r.auth_seq_id is not None and r.ca_coords is not None
        }

        seq_pos = {seq: pos for pos, seq in enumerate(residues)}

        # Build adjacency: contact if distance < threshold, excluding
        # sequence-adjacent pairs (|i - j| <= 1).
        adjacency: Dict[int, List[int]] = {seq: [] for seq in residues}
        seqs = residues
        for a_pos in range(len(seqs)):
            for b_pos in range(a_pos + 1, len(seqs)):
                seq_a, seq_b = seqs[a_pos], seqs[b_pos]
                if abs(seq_pos[seq_a] - seq_pos[seq_b]) <= 1:
                    continue
                dist = float(np.linalg.norm(coords[seq_a] - coords[seq_b]))
                if dist < contact_threshold:
                    adjacency[seq_a].append(seq_b)
                    adjacency[seq_b].append(seq_a)

        # Connected components (BFS).
        visited = set()
        components: List[List[int]] = []
        for seq in seqs:
            if seq in visited:
                continue
            stack = [seq]
            component = []
            visited.add(seq)
            while stack:
                current = stack.pop()
                component.append(current)
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        stack.append(neighbor)
            components.append(sorted(component))

        for component in components:
            if len(component) < min_residues:
                continue
            regions.append({
                "chain": chain_id,
                "start": min(component),
                "end": max(component),
                "seq_list": list(component),
                "n_residues": len(component),
            })

    # Stable ordering: chain, then start; labels assigned last.
    regions.sort(key=lambda r: (r["chain"], r["start"]))
    for k, region in enumerate(regions, start=1):
        region["label"] = f"region_{k}"

    return regions


def _robust_z_scores(scores: np.ndarray) -> Optional[np.ndarray]:
    """Robust z-scores for flagging displacement peaks.

    Primary: z = (x - median) / (1.4826 * MAD).

    Degenerate case (MAD == 0, i.e. more than half the residues share one
    identical score — the 'bulk' level): every residue whose score differs
    from the bulk level is maximally anomalous relative to a zero-spread
    bulk and is returned with an infinite z. Returns None only when all
    scores are identical (no spread at all, nothing to flag).
    """
    med = float(np.median(scores))
    mad = float(np.median(np.abs(scores - med)))
    if mad > 0:
        return (scores - med) / (1.4826 * mad)
    if not np.any(scores != med):
        return None
    z = np.zeros_like(scores)
    z[scores != med] = np.inf
    return z


def _robust_outlier_indices(
    scores: np.ndarray,
    z_threshold: float,
) -> List[int]:
    """Indices whose robust z-score exceeds the threshold."""
    z = _robust_z_scores(scores)
    if z is None:
        return []
    return [int(i) for i in np.where(z >= z_threshold)[0]]


__all__ = [
    "detect_displacement_sites",
    "detect_contact_regions",
]
