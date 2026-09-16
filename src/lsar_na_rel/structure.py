"""8ACT biological-assembly structure features (Talk-v0 / #28 Task 3).

Model columns (only):

- ``wt_site_rsa_8act_mean_ab`` — arithmetic mean of chain-A/B DSSP RSA
  (ASA / Sander maximum via Biopython ``acc_array="Sander"``)
- ``ihm_interface_flag_8act`` — OR of BH–FH and BH–S2 geometry subtype bits
  at ≤4.5 Å non-hydrogen contacts

Chain roles: BH=A, FH=B (Grinzato R403(A)–Y455(B) contact). Proximal S2 is a
residue-range resegmentation on A/B (no separate S2 chain). Light chains C–F
do not set the IHM flag. Volume/hydropathy/contact-count/FoldX columns are
not generated. Feature generation accepts variant labels only.
"""

from __future__ import annotations

import gzip
import hashlib
import os
import tempfile
import urllib.request
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping, Sequence

import numpy as np
import pandas as pd
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import DSSP, PDBParser
from Bio.PDB.Polypeptide import is_aa

from lsar_na_rel.inputs import parse_variant

ASSEMBLY_URL = "https://files.rcsb.org/download/8ACT.pdb1.gz"
GZIP_SHA256 = "84469940f04dd8428a98e76c553ab519c1d45aef2beda176b8813d3b120e797b"
DECOMPRESSED_PDB_SHA256 = (
    "ef548f824d10986e313fa41910fa10d5d3925bc04ae3350f58214a2639d3d6a9"
)

CONTACT_CUTOFF_ANGSTROM = 4.5
DOWNLOAD_TIMEOUT_SECONDS = 60.0
BH_CHAIN = "A"
FH_CHAIN = "B"
HEAVY_CHAINS = (BH_CHAIN, FH_CHAIN)
LIGHT_CHAINS = ("C", "D", "E", "F")
REQUIRED_CHAINS = frozenset({"A", "B", "C", "D", "E", "F"})

# Frozen S1/S2 sequence cut (symmetric on A/B), bounded by modeled 8ACT residues.
# Primary sources: PDB 2FXM COMPND/DBREF (human cardiac β-myosin S2-delta =
# P12883 838–963); Blankenfeldt et al. state invariant Pro838 marks the
# C terminus of S1. Grinzato 8ACT is MYH7 1–942 with 15 heptads of proximal S2.
# Deposited continuous coiled-coil helix onsets A841/B846 are asymmetric
# coiled-coil formation, not the S1/S2 sequence boundary — retained as
# provenance context only.
MOTOR_LEVER_RANGES: dict[str, tuple[int, int]] = {
    "A": (3, 837),
    "B": (3, 837),
}
PROXIMAL_S2_RANGES: dict[str, tuple[int, int]] = {
    "A": (838, 906),
    "B": (838, 898),
}
DEPOSITED_HELIX_ONSETS: dict[str, int] = {
    "A": 841,
    "B": 846,
}
DEPOSITED_HELIX_RANGES: dict[str, tuple[int, int]] = {
    "A": (841, 906),
    "B": (846, 898),
}

# Modified lysine used only for WT identity normalization.
_M3L_TO_CANONICAL = "K"
_NONSTANDARD_IDENTITY = {
    "M3L": _M3L_TO_CANONICAL,
}


@dataclass(frozen=True)
class AssemblyManifest:
    """Checksummed local 8ACT biological-assembly cache record."""

    url: str
    cache_dir: Path
    gzip_path: Path
    pdb_path: Path
    gzip_sha256: str
    pdb_sha256: str


@dataclass(frozen=True)
class StructureResult:
    """Model feature table plus row-level and assembly provenance."""

    features: pd.DataFrame
    site_provenance: pd.DataFrame
    provenance: Mapping[str, object]


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _download_url_to_path(url: str, target: Path) -> None:
    """Download ``url`` into ``target`` via a temporary partial file."""
    target.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{target.name}.",
        suffix=".partial",
        dir=str(target.parent),
    )
    tmp_path = Path(tmp_name)
    try:
        with os.fdopen(fd, "wb") as handle, urllib.request.urlopen(
            url,
            timeout=DOWNLOAD_TIMEOUT_SECONDS,
        ) as response:
            while True:
                chunk = response.read(1024 * 1024)
                if not chunk:
                    break
                handle.write(chunk)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(tmp_path, target)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        # Never leave an unverified final artifact after a failed download.
        target.unlink(missing_ok=True)
        raise


def _decompress_gzip_atomic(gzip_path: Path, pdb_path: Path) -> None:
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{pdb_path.name}.",
        suffix=".partial",
        dir=str(pdb_path.parent),
    )
    tmp_path = Path(tmp_name)
    try:
        with gzip.open(gzip_path, "rb") as src, os.fdopen(fd, "wb") as dst:
            while True:
                chunk = src.read(1024 * 1024)
                if not chunk:
                    break
                dst.write(chunk)
            dst.flush()
            os.fsync(dst.fileno())
        os.replace(tmp_path, pdb_path)
    except Exception:
        tmp_path.unlink(missing_ok=True)
        raise


def ensure_8act_assembly(
    cache_dir: str | Path,
    *,
    url: str = ASSEMBLY_URL,
    gzip_sha256: str = GZIP_SHA256,
    pdb_sha256: str = DECOMPRESSED_PDB_SHA256,
    gzip_name: str = "8ACT.pdb1.gz",
    pdb_name: str = "8ACT-assembly1.pdb",
) -> AssemblyManifest:
    """Ensure a checksum-matched biological assembly is present locally.

    Downloads ``url`` when the gzip is absent, verifies gzip bytes before
    decompression, then verifies the decompressed PDB. Uses atomic writes and
    removes partial files on failure.
    """
    cache = Path(cache_dir)
    cache.mkdir(parents=True, exist_ok=True)
    gzip_path = cache / gzip_name
    pdb_path = cache / pdb_name

    try:
        if not gzip_path.is_file():
            _download_url_to_path(url, gzip_path)

        observed_gzip = _sha256_file(gzip_path)
        if observed_gzip != gzip_sha256:
            gzip_path.unlink(missing_ok=True)
            pdb_path.unlink(missing_ok=True)
            raise ValueError(
                f"Assembly gzip SHA-256 mismatch for {gzip_path}: "
                f"{observed_gzip} != {gzip_sha256}"
            )

        if not pdb_path.is_file():
            _decompress_gzip_atomic(gzip_path, pdb_path)

        observed_pdb = _sha256_file(pdb_path)
        if observed_pdb != pdb_sha256:
            # Repair from the already-verified gzip in the same call.
            pdb_path.unlink(missing_ok=True)
            _decompress_gzip_atomic(gzip_path, pdb_path)
            observed_pdb = _sha256_file(pdb_path)
            if observed_pdb != pdb_sha256:
                pdb_path.unlink(missing_ok=True)
                raise ValueError(
                    f"Assembly PDB SHA-256 mismatch for {pdb_path}: "
                    f"{observed_pdb} != {pdb_sha256}"
                )
    except Exception:
        # Best-effort cleanup of temp partials left in the cache directory.
        for leftover in cache.glob("*.partial"):
            leftover.unlink(missing_ok=True)
        for leftover in cache.glob("*.tmp"):
            leftover.unlink(missing_ok=True)
        raise

    return AssemblyManifest(
        url=url,
        cache_dir=cache,
        gzip_path=gzip_path,
        pdb_path=pdb_path,
        gzip_sha256=observed_gzip,
        pdb_sha256=observed_pdb,
    )


def verify_assembly_pdb_sha256(
    pdb_path: str | Path,
    *,
    expected_sha256: str = DECOMPRESSED_PDB_SHA256,
) -> str:
    """Hash ``pdb_path`` and require it match the pinned decompressed digest."""
    path = Path(pdb_path)
    if not path.is_file():
        raise FileNotFoundError(f"Assembly PDB not found: {path}")
    observed = _sha256_file(path)
    if observed != expected_sha256:
        raise ValueError(
            f"Assembly PDB SHA-256 mismatch for {path}: "
            f"{observed} != {expected_sha256}"
        )
    return observed


def load_8act_model(pdb_path: str | Path):
    """Parse the full biological assembly and require chains A–F (any order)."""
    path = Path(pdb_path)
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("8ACT", str(path))
    model = structure[0]
    chain_ids = {chain.id for chain in model}
    if chain_ids != set(REQUIRED_CHAINS):
        raise ValueError(
            f"Expected chains {sorted(REQUIRED_CHAINS)}; found {sorted(chain_ids)}"
        )
    return model


def normalize_residue_identity(resname: str) -> str:
    """Map a PDB residue name to a one-letter code (M3L → K)."""
    code = str(resname).strip().upper()
    if code in _NONSTANDARD_IDENTITY:
        return _NONSTANDARD_IDENTITY[code]
    if len(code) == 1:
        return code
    key = code[0] + code[1:].lower()
    try:
        return protein_letters_3to1[key].upper()
    except KeyError as exc:
        raise ValueError(f"Unrecognized residue name for identity: {resname!r}") from exc


def _residue_at(model, chain_id: str, position: int):
    chain = model[chain_id]
    # Prefer standard polymeric residue at this sequence number.
    key = (" ", position, " ")
    if key in chain:
        residue = chain[key]
        if is_aa(residue, standard=False):
            return residue
    matches = [
        residue
        for residue in chain
        if residue.id[1] == position and is_aa(residue, standard=False)
    ]
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(f"Ambiguous residues at {chain_id}:{position}")
    return matches[0]


def validate_wt_sites_ab(model, variant_labels: Sequence[str]) -> dict[str, dict[str, object]]:
    """Validate claimed WT identity on heavy chains A and B for each variant."""
    sites: dict[str, dict[str, object]] = {}
    for label in variant_labels:
        parsed = parse_variant(str(label))
        observed: dict[str, object] = {
            "claimed_wt": parsed.wt,
            "position": parsed.position,
        }
        for chain_id in HEAVY_CHAINS:
            residue = _residue_at(model, chain_id, parsed.position)
            if residue is None:
                raise ValueError(
                    f"Missing WT residue for {label!r} on chain {chain_id} "
                    f"at position {parsed.position}"
                )
            identity = normalize_residue_identity(residue.get_resname())
            if identity != parsed.wt:
                raise ValueError(
                    f"WT identity mismatch for {label!r} on chain {chain_id}: "
                    f"claimed {parsed.wt!r} but assembly has {identity!r} "
                    f"({residue.get_resname().strip()})"
                )
            observed[chain_id] = identity
            observed[f"covered_{chain_id.lower()}"] = True
        sites[parsed.label] = observed
    return sites


def heavy_atom_contact(
    coord_a: Sequence[float],
    coord_b: Sequence[float],
    *,
    cutoff: float = CONTACT_CUTOFF_ANGSTROM,
) -> bool:
    """True when Euclidean distance between two points is ≤ cutoff."""
    a = np.asarray(coord_a, dtype=float)
    b = np.asarray(coord_b, dtype=float)
    return float(np.linalg.norm(a - b)) <= float(cutoff)


def _is_hydrogen_like_name(name: str) -> bool:
    """True for PDB hydrogen/deuterium names, including digit-prefixed forms."""
    token = name.strip().upper()
    if not token:
        return False
    # Strip leading digits used in PDB hydrogen naming (e.g. 1HB, 2HD1).
    i = 0
    while i < len(token) and token[i].isdigit():
        i += 1
    core = token[i:]
    if not core:
        return False
    return core[0] in {"H", "D"}


def _is_hydrogen_atom(atom) -> bool:
    element = (atom.element or "").strip().upper()
    if element in {"H", "D"}:
        return True
    # Guard PDB names when element is missing/blank (incl. digit-prefixed H/D).
    if element in {"", " "}:
        return _is_hydrogen_like_name(atom.get_name())
    return False


def residues_in_heavy_contact(residue_a, residue_b, *, cutoff: float = CONTACT_CUTOFF_ANGSTROM) -> bool:
    """True if any non-hydrogen atom pair is within ``cutoff`` Å."""
    atoms_a = [atom for atom in residue_a.get_atoms() if not _is_hydrogen_atom(atom)]
    atoms_b = [atom for atom in residue_b.get_atoms() if not _is_hydrogen_atom(atom)]
    cutoff_sq = float(cutoff) ** 2
    for atom_a in atoms_a:
        coord_a = atom_a.coord
        for atom_b in atoms_b:
            delta = coord_a - atom_b.coord
            if float(np.dot(delta, delta)) <= cutoff_sq:
                return True
    return False


def _in_inclusive_range(position: int, span: tuple[int, int]) -> bool:
    start, end = span
    return start <= position <= end


def ihm_partner_contacts(
    model,
    position: int,
    *,
    cutoff: float = CONTACT_CUTOFF_ANGSTROM,
) -> dict[str, int]:
    """Geometry-only BH–FH / BH–S2 flags for the WT site on A and/or B."""
    bh_fh = 0
    bh_s2 = 0

    for focal_chain in HEAVY_CHAINS:
        if focal_chain not in model:
            continue
        focal = _residue_at(model, focal_chain, position)
        if focal is None:
            continue
        partner_chain = FH_CHAIN if focal_chain == BH_CHAIN else BH_CHAIN
        if partner_chain not in model:
            continue
        partner = model[partner_chain]
        motor = MOTOR_LEVER_RANGES[partner_chain]
        s2 = PROXIMAL_S2_RANGES[partner_chain]
        for residue in partner:
            if not is_aa(residue, standard=False):
                continue
            resseq = int(residue.id[1])
            if not residues_in_heavy_contact(focal, residue, cutoff=cutoff):
                continue
            if _in_inclusive_range(resseq, motor):
                bh_fh = 1
            elif _in_inclusive_range(resseq, s2):
                bh_s2 = 1
        # Light-chain contacts are intentionally ignored for both subtypes.

    return {
        "ihm_contact_bh_fh": bh_fh,
        "ihm_contact_bh_s2": bh_s2,
        "ihm_interface_flag_8act": 1 if (bh_fh or bh_s2) else 0,
    }


def _rsa_from_dssp(dssp, positions: Sequence[int]) -> dict[tuple[str, int], float]:
    """Extract Sander-relative accessibility for requested A/B positions."""
    wanted = {(chain, int(pos)) for chain in HEAVY_CHAINS for pos in positions}
    values: dict[tuple[str, int], float] = {}
    for key in dssp.keys():
        chain_id, res_id = key[0], key[1]
        if chain_id not in HEAVY_CHAINS:
            continue
        resseq = int(res_id[1])
        pair = (chain_id, resseq)
        if pair not in wanted:
            continue
        # Biopython DSSP tuple: index 3 is relative ASA under acc_array.
        values[pair] = float(dssp[key][3])
    missing = wanted - set(values)
    if missing:
        raise ValueError(f"DSSP missing RSA for sites: {sorted(missing)}")
    return values


def generate_structure_features(
    variant_labels: Sequence[str],
    *,
    assembly_pdb_path: str | Path,
    mkdssp: str = "mkdssp",
    contact_cutoff: float = CONTACT_CUTOFF_ANGSTROM,
    expected_pdb_sha256: str = DECOMPRESSED_PDB_SHA256,
) -> StructureResult:
    """Compute Talk-v0 structure features from variant labels (no outcomes)."""
    pdb_path = Path(assembly_pdb_path)
    # Provenance honesty: verify the bytes actually used before DSSP/features.
    observed_pdb_sha256 = verify_assembly_pdb_sha256(
        pdb_path,
        expected_sha256=expected_pdb_sha256,
    )
    model = load_8act_model(pdb_path)
    sites = validate_wt_sites_ab(model, variant_labels)
    positions = [int(info["position"]) for info in sites.values()]

    with warnings.catch_warnings():
        # mkdssp may emit a harmless PDB/mmCIF probe warning on PDB input.
        warnings.filterwarnings(
            "ignore",
            message=r".*does not seem to be an mmCIF file.*",
            category=UserWarning,
        )
        dssp = DSSP(model, str(pdb_path), dssp=mkdssp, acc_array="Sander")
    rsa = _rsa_from_dssp(dssp, positions)

    feature_rows: list[dict[str, object]] = []
    provenance_rows: list[dict[str, object]] = []
    for label, info in sites.items():
        position = int(info["position"])
        rsa_a = float(rsa[(BH_CHAIN, position)])
        rsa_b = float(rsa[(FH_CHAIN, position)])
        mean_ab = (rsa_a + rsa_b) / 2.0
        flags = ihm_partner_contacts(model, position, cutoff=contact_cutoff)
        feature_rows.append(
            {
                "variant": label,
                "wt_site_rsa_8act_mean_ab": mean_ab,
                "ihm_interface_flag_8act": int(flags["ihm_interface_flag_8act"]),
            }
        )
        provenance_rows.append(
            {
                "variant": label,
                "position": position,
                "wt_site_rsa_8act_a": rsa_a,
                "wt_site_rsa_8act_b": rsa_b,
                "coverage_a": True,
                "coverage_b": True,
                "wt_site_rsa_8act_bh": rsa_a,
                "wt_site_rsa_8act_fh": rsa_b,
                "ihm_contact_bh_fh": int(flags["ihm_contact_bh_fh"]),
                "ihm_contact_bh_s2": int(flags["ihm_contact_bh_s2"]),
            }
        )

    features = pd.DataFrame(
        feature_rows,
        columns=["variant", "wt_site_rsa_8act_mean_ab", "ihm_interface_flag_8act"],
    )
    site_provenance = pd.DataFrame(provenance_rows)
    provenance: dict[str, object] = {
        "assembly_url": ASSEMBLY_URL,
        "assembly_pdb_path": str(pdb_path),
        # Observed/verified for this generation path (PDB bytes only).
        "pdb_sha256": observed_pdb_sha256,
        "pdb_sha256_verified": True,
        "gzip_sha256_verified": False,
        # Source/pin metadata (not claimed as verified on this path).
        "pinned_pdb_sha256": DECOMPRESSED_PDB_SHA256,
        "pinned_gzip_sha256": GZIP_SHA256,
        "bh_chain": BH_CHAIN,
        "fh_chain": FH_CHAIN,
        "motor_lever_ranges": dict(MOTOR_LEVER_RANGES),
        "proximal_s2_ranges": dict(PROXIMAL_S2_RANGES),
        "deposited_helix_onsets": dict(DEPOSITED_HELIX_ONSETS),
        "deposited_helix_ranges": dict(DEPOSITED_HELIX_RANGES),
        "s1_s2_cut_rationale": (
            "Symmetric sequence cut motor/lever<=837, proximal S2>=838 "
            "(PDB 2FXM / Blankenfeldt Pro838 S1 C-terminus; Grinzato 8ACT "
            "MYH7 1-942 proximal S2). Deposited helix onsets A841/B846 are "
            "coiled-coil formation provenance, not classifier boundaries."
        ),
        "contact_cutoff_angstrom": float(contact_cutoff),
        "dssp": mkdssp,
        "acc_array": "Sander",
        "feature_columns": (
            "wt_site_rsa_8act_mean_ab",
            "ihm_interface_flag_8act",
        ),
        "excluded_model_columns": (
            "delta_volume",
            "delta_hydropathy",
            "contact_count_8act_mean_ab",
            "foldx_ddg_bh_fh",
            "secondary_structure",
            "domain",
        ),
        "rsa_definition": "(RSA_A + RSA_B) / 2 with RSA = DSSP ASA / Sander maximum",
        "ihm_definition": (
            "OR of ihm_contact_bh_fh and ihm_contact_bh_s2; "
            "≤4.5 Å non-hydrogen contacts; light chains C–F excluded"
        ),
    }
    return StructureResult(
        features=features,
        site_provenance=site_provenance,
        provenance=provenance,
    )
