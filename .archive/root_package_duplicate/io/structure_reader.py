"""
mmCIF structure reader for AF3 structural analysis.

Parses AF3-predicted mmCIF files into NormalisedStructure objects.
Uses gemmi for mmCIF parsing; no biological assumptions are made.

The rest of the pipeline operates on NormalisedStructure, not mmCIF.
"""

from __future__ import annotations

import hashlib
from pathlib import Path
from typing import List, Optional, Tuple

from af3_analysis.io.provenance import Checksum
from af3_analysis.structural.enums import ParseStatus
from af3_analysis.structural.representation import (
    AtomRecord,
    ChainInfo,
    EntityInfo,
    NormalisedStructure,
    ResidueRecord,
)


class StructureParseError(Exception):
    """Error parsing an mmCIF structure file."""
    pass


def _compute_checksum(path: Path) -> str:
    """Compute SHA-256 checksum of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _safe_int(value: str) -> Optional[int]:
    """Parse an integer from mmCIF, returning None for '.' or '?'."""
    value = value.strip()
    if value in (".", "?", ""):
        return None
    try:
        return int(value)
    except ValueError:
        return None


def _safe_float(value: str) -> Optional[float]:
    """Parse a float from mmCIF, returning None for '.' or '?'."""
    value = value.strip()
    if value in (".", "?", ""):
        return None
    try:
        return float(value)
    except ValueError:
        return None


def parse_mmcif(
    cif_path: Path,
    condition_id: str,
    seed: int,
    sample: int,
    *,
    model_num: Optional[int] = None,
) -> NormalisedStructure:
    """
    Parse an mmCIF file into a NormalisedStructure.

    Parameters
    ----------
    cif_path : Path
        Path to the mmCIF file.
    condition_id : str
        Condition identifier for this prediction.
    seed : int
        AF3 seed value.
    sample : int
        AF3 sample value.
    model_num : int, optional
        Model number to extract. If None, uses the first model found.

    Returns
    -------
    NormalisedStructure
        Parsed structure with atoms, residues, entities, chains.

    Raises
    ------
    StructureParseError
        If the file cannot be parsed or contains no atoms.
    """
    import gemmi

    try:
        doc = gemmi.cif.read_file(str(cif_path))
    except Exception as e:
        raise StructureParseError(f"Cannot parse mmCIF file {cif_path}: {e}")

    block = doc.sole_block()

    # --- Extract entity metadata ---
    entities = _extract_entities(block)
    chains = _extract_chains(block, entities)

    # --- Extract atom_site data ---
    atoms = _extract_atoms(block, model_num)

    if not atoms:
        raise StructureParseError(
            f"No atoms found in {cif_path} (model_num={model_num})"
        )

    # --- Build residue records ---
    residues = _build_residues(atoms)

    # --- Compute checksum ---
    checksum = _compute_checksum(cif_path)

    # --- Build NormalisedStructure ---
    structure = NormalisedStructure(
        source_path=cif_path,
        source_checksum=checksum,
        condition_id=condition_id,
        seed=seed,
        sample=sample,
        atoms=atoms,
        residues=residues,
        entities=entities,
        chains=chains,
    )
    structure.build_indices()

    return structure


def _extract_entities(block) -> List[EntityInfo]:
    """Extract entity metadata from mmCIF block."""
    import gemmi

    entities = []

    # Get entity info using find_values (more robust across gemmi versions)
    entity_ids = []
    entity_types = []
    entity_descriptions = []

    if block.find_values("_entity.id"):
        entity_ids = [str(v).strip() for v in block.find_values("_entity.id")]
    if block.find_values("_entity.type"):
        entity_types = [str(v).strip() for v in block.find_values("_entity.type")]
    if block.find_values("_entity.pdbx_description"):
        entity_descriptions = [str(v).strip() for v in block.find_values("_entity.pdbx_description")]

    # Get polymer types using find_values
    polymer_types = {}  # entity_id -> polymer_type
    if block.find_values("_entity_poly.entity_id"):
        ep_eids = [str(v).strip() for v in block.find_values("_entity_poly.entity_id")]
        ep_types = [str(v).strip() for v in block.find_values("_entity_poly.type")]
        for j, eid in enumerate(ep_eids):
            ptype = ep_types[j] if j < len(ep_types) else None
            polymer_types[eid] = ptype

    # Get chain-to-entity mapping from _struct_asym using find_values
    chain_entities = {}  # chain_id -> entity_id
    chain_auth = {}  # chain_id -> auth_asym_id
    if block.find_values("_struct_asym.entity_id"):
        sa_eids = [str(v).strip() for v in block.find_values("_struct_asym.entity_id")]
        sa_ids = [str(v).strip() for v in block.find_values("_struct_asym.id")]
        for j, eid in enumerate(sa_eids):
            cid = sa_ids[j] if j < len(sa_ids) else ""
            chain_entities[cid] = eid
            chain_auth[cid] = cid  # default: label == auth

    # Build entity info
    for i, eid in enumerate(entity_ids):
        etype = entity_types[i] if i < len(entity_types) else "unknown"
        desc = entity_descriptions[i] if i < len(entity_descriptions) else None
        ptype = polymer_types.get(eid)

        # Find the chain for this entity
        chain_id = None
        for cid, ceid in chain_entities.items():
            if ceid == eid:
                chain_id = cid
                break

        entities.append(EntityInfo(
            entity_id=int(eid),
            chain_id=chain_id or "",
            auth_asym_id=chain_auth.get(chain_id, chain_id or ""),
            entity_type=etype,
            polymer_type=ptype,
            description=desc,
            n_residues=0,
            n_atoms=0,
        ))

    return entities


def _extract_chains(block, entities: List[EntityInfo]) -> List[ChainInfo]:
    """Extract chain metadata from mmCIF block."""
    chains = []
    for entity in entities:
        chains.append(ChainInfo(
            chain_id=entity.chain_id,
            auth_asym_id=entity.auth_asym_id,
            entity_id=entity.entity_id,
            entity_type=entity.entity_type,
            polymer_type=entity.polymer_type,
            n_residues=entity.n_residues,
            n_atoms=entity.n_atoms,
            first_seq_id=None,
            last_seq_id=None,
        ))
    return chains


def _extract_atoms(block, model_num: Optional[int] = None) -> List[AtomRecord]:
    """Extract atom coordinates from mmCIF _atom_site loop."""
    atoms = []

    if not block.find_values("_atom_site.id"):
        return atoms

    # Determine model number to use
    if model_num is None:
        # Use the first model number found
        for row in block.find_loop("_atom_site.group_PDB"):
            # row contains all columns for this atom
            # We need the full row, so use _atom_site.id to get all rows
            break
        # Actually, find_loop with a single column name returns one column.
        # We need to iterate using find_values for the full loop.
        # Use _atom_site.id to discover the number of atoms, then parse.
        id_values = list(block.find_values("_atom_site.id"))
        if id_values:
            model_num = 1  # Default; refined below if pdbx_PDB_model_num column exists

    # Parse all atom rows using find_values for each column
    # This is more robust than find_loop across gemmi versions
    group_pdb_vals = list(block.find_values("_atom_site.group_PDB"))
    id_vals = list(block.find_values("_atom_site.id"))
    type_vals = list(block.find_values("_atom_site.type_symbol"))
    atom_name_vals = list(block.find_values("_atom_site.label_atom_id"))
    comp_id_vals = list(block.find_values("_atom_site.label_comp_id"))
    asym_vals = list(block.find_values("_atom_site.label_asym_id"))
    entity_vals = list(block.find_values("_atom_site.label_entity_id"))
    seq_vals = list(block.find_values("_atom_site.label_seq_id"))
    x_vals = list(block.find_values("_atom_site.Cartn_x"))
    y_vals = list(block.find_values("_atom_site.Cartn_y"))
    z_vals = list(block.find_values("_atom_site.Cartn_z"))
    occ_vals = list(block.find_values("_atom_site.occupancy"))
    b_vals = list(block.find_values("_atom_site.B_iso_or_equiv"))
    auth_seq_vals = list(block.find_values("_atom_site.auth_seq_id"))
    auth_asym_vals = list(block.find_values("_atom_site.auth_asym_id"))
    model_vals = list(block.find_values("_atom_site.pdbx_PDB_model_num"))

    n_atoms = len(group_pdb_vals)
    if n_atoms == 0:
        return atoms

    # Determine model number from first row
    if model_num is None and model_vals:
        model_num = _safe_int(str(model_vals[0]).strip()) or 1
    elif model_num is None:
        model_num = 1

    # Parse all atom rows
    for i in range(n_atoms):
        group_pdb = str(group_pdb_vals[i]).strip()

        # Only ATOM and HETATM records
        if group_pdb not in ("ATOM", "HETATM"):
            continue

        # Parse fields from column arrays
        atom_id = _safe_int(str(id_vals[i]).strip()) if i < len(id_vals) else None
        atom_type = str(type_vals[i]).strip() if i < len(type_vals) else ""
        atom_name = str(atom_name_vals[i]).strip() if i < len(atom_name_vals) else ""
        comp_id = str(comp_id_vals[i]).strip() if i < len(comp_id_vals) else ""
        chain_id = str(asym_vals[i]).strip() if i < len(asym_vals) else ""
        entity_id = _safe_int(str(entity_vals[i]).strip()) if i < len(entity_vals) else 1
        seq_id = _safe_int(str(seq_vals[i]).strip()) if i < len(seq_vals) else None
        x = _safe_float(str(x_vals[i]).strip()) if i < len(x_vals) else None
        y = _safe_float(str(y_vals[i]).strip()) if i < len(y_vals) else None
        z = _safe_float(str(z_vals[i]).strip()) if i < len(z_vals) else None
        occupancy = _safe_float(str(occ_vals[i]).strip()) if i < len(occ_vals) else 1.0
        b_factor = _safe_float(str(b_vals[i]).strip()) if i < len(b_vals) else None
        auth_seq_id = _safe_int(str(auth_seq_vals[i]).strip()) if i < len(auth_seq_vals) else None
        auth_asym_id = str(auth_asym_vals[i]).strip() if i < len(auth_asym_vals) else chain_id
        row_model_num = _safe_int(str(model_vals[i]).strip()) if i < len(model_vals) else 1

        # Skip atoms not in the requested model
        if row_model_num is not None and row_model_num != model_num:
            continue

        # Skip atoms without coordinates
        if x is None or y is None or z is None:
            continue

        if atom_id is None:
            continue

        atoms.append(AtomRecord(
            atom_id=atom_id,
            atom_name=atom_name,
            atom_type=atom_type,
            comp_id=comp_id,
            chain_id=chain_id,
            entity_id=entity_id or 1,
            seq_id=seq_id,
            auth_seq_id=auth_seq_id,
            auth_asym_id=auth_asym_id,
            coords=(x, y, z),
            b_factor=b_factor,
            occupancy=occupancy,
            model_num=row_model_num or model_num,
        ))

    return atoms


def _build_residues(atoms: List[AtomRecord]) -> List[ResidueRecord]:
    """Build residue-level records from atom records."""
    from collections import defaultdict

    residue_atoms = defaultdict(list)
    for atom in atoms:
        key = (atom.chain_id, atom.auth_seq_id or atom.seq_id)
        residue_atoms[key].append(atom)

    residues = []
    for (chain_id, seq_id), res_atoms in residue_atoms.items():
        atom_names = tuple(sorted(set(a.atom_name for a in res_atoms)))
        n_atoms = len(res_atoms)
        comp_id = res_atoms[0].comp_id
        entity_id = res_atoms[0].entity_id

        # Check completeness (protein backbone)
        is_complete = False
        if res_atoms[0].entity_id:
            # Heuristic: if we have N, CA, C atoms, it's likely complete
            name_set = set(atom_names)
            if {"N", "CA", "C"}.issubset(name_set):
                is_complete = True

        residues.append(ResidueRecord(
            comp_id=comp_id,
            chain_id=chain_id,
            entity_id=entity_id,
            seq_id=None,
            auth_seq_id=seq_id,
            atom_names=atom_names,
            n_atoms=n_atoms,
            is_complete=is_complete,
        ))

    # Sort by chain and sequence
    residues.sort(key=lambda r: (r.chain_id, r.auth_seq_id or 0))

    return residues


def discover_structure_files(
    root: Path,
    *,
    pattern: str = "*_model.cif",
) -> List[Tuple[Path, str, int, int]]:
    """
    Discover AF3 structure files under a root directory.

    Returns list of (path, condition_id, seed, sample) tuples.
    Uses the same naming convention as RawAF3Reader.MODEL_PATTERNS.
    """
    import re

    model_pattern = re.compile(
        r"^(?P<cond>.+)_seed-(?P<seed>\d+)_sample-(?P<sample>\d+)_model\.cif$"
    )

    results = []
    for cif_path in root.rglob(pattern):
        if not cif_path.is_file():
            continue

        match = model_pattern.match(cif_path.name)
        if match:
            cond = match.group("cond")
            seed = int(match.group("seed"))
            sample = int(match.group("sample"))
            results.append((cif_path, cond, seed, sample))

    return sorted(results, key=lambda x: (x[1], x[2], x[3]))
