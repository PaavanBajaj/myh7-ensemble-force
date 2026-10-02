"""Label-blind, checksum-pinned WT geometry for the MYH7 velocity screen."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from Bio import __version__ as biopython_version
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict


EXPECTED_SHA256 = {
    "fasta": "45c96586dd50e0ede770a3bb4fb0aab125d103d53b2bc0b0fed89619e07a16ad",
    "8ACT": "ef548f824d10986e313fa41910fa10d5d3925bc04ae3350f58214a2639d3d6a9",
    "8EFE": "12d0d9ac570f2775fe05144d09a2e15763df5daf1478d6a9b5d3623e224a4800",
    "8EFD": "fed8cd752b452399b8409ee72dfa8785625f3a6cd9fb5919f8e961aa94e12e58",
    "8EFI": "ffd2c42aeef49d1d27ad8e6a71826eca9b2ad95b119d375d8a5acb052ab9dfe8",
}
STRUCTURE_URLS = {
    "8ACT": "https://files.rcsb.org/download/8ACT.pdb1.gz",
    "8EFE": "https://files.rcsb.org/download/8EFE.cif",
    "8EFD": "https://files.rcsb.org/download/8EFD.cif",
    "8EFI": "https://files.rcsb.org/download/8EFI.cif",
}
FEATURES = (
    "adp_rigor_ca_shift_A",
    "adp_adenine_ca_distance_8efe_A",
    "rigor_actin_min_heavy_distance_8efi_A",
    "relay_landmark_min_ca_distance_8act_mean_ab",
    "elc_min_heavy_distance_8act_mean_ab",
)
RELAY_LANDMARKS = (497, 501, 504, 506, 508, 510, 511, 513)
ALIGNMENT_CORE = ((6, 199), (216, 624))
ADENINE_ATOMS = frozenset(("N9", "C8", "N7", "C5", "C6", "N6", "N1", "C2", "N3", "C4"))
VARIANT_RE = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])([1-9][0-9]*)([ACDEFGHIKLMNPQRSTVWY])$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_input(name: str, path: Path) -> str:
    observed = sha256(path)
    if observed != EXPECTED_SHA256[name]:
        raise ValueError(f"{name} SHA-256 mismatch: {observed} != {EXPECTED_SHA256[name]}")
    return observed


def residue_letter(name: str) -> str:
    if name.upper() == "M3L":
        return "K"
    try:
        return protein_letters_3to1[name.title()].upper()
    except KeyError as exc:
        raise ValueError(f"unrecognized amino acid {name!r}") from exc


def load_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError("expected FASTA with a header")
    sequence = "".join(line.strip() for line in lines[1:]).upper()
    if len(sequence) != 1935 or hashlib.sha256(sequence.encode()).hexdigest() != "e98e5d01a359820bcd1a00c3422025d9fedf1fdbc46ab586e3f5b731f4585be6":
        raise ValueError("P12883 sequence identity/checksum mismatch")
    return sequence


def load_variants(path: Path, sequence: str) -> list[tuple[str, int, str]]:
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or set(rows[0]) != {"variant"}:
        raise ValueError("variant registry must have only a variant column")
    result = []
    seen = set()
    for row in rows:
        label = row["variant"]
        match = VARIANT_RE.fullmatch(label)
        if match is None or label in seen:
            raise ValueError(f"invalid or duplicate variant {label!r}")
        wt, number, mutant = match.groups()
        position = int(number)
        if position > len(sequence) or sequence[position - 1] != wt or wt == mutant:
            raise ValueError(f"variant {label} does not match P12883")
        result.append((label, position, wt))
        seen.add(label)
    return result


@dataclass
class Site:
    residue: str
    auth_chain: str
    auth_seq: int
    label_seq: int
    ca: np.ndarray | None
    atoms: np.ndarray


@dataclass
class CifStructure:
    title: str
    sites: dict[str, dict[int, Site]]
    ligand_atoms: dict[str, np.ndarray]


def _coordinate(columns: dict[str, list[str]], i: int) -> np.ndarray:
    return np.array([float(columns[f"_atom_site.Cartn_{axis}"][i]) for axis in "xyz"], dtype=float)


def load_cif(path: Path, sequence: str, *, chains: set[str]) -> CifStructure:
    data = MMCIF2Dict(str(path))
    if "P12883" not in data.get("_struct_ref.pdbx_db_accession", []):
        raise ValueError(f"{path.name} lacks P12883 polymer reference")
    fields = (
        "group_PDB", "label_atom_id", "label_comp_id", "label_asym_id", "auth_asym_id",
        "label_seq_id", "auth_seq_id", "type_symbol", "occupancy", "label_alt_id",
        "pdbx_PDB_model_num", "Cartn_x", "Cartn_y", "Cartn_z",
    )
    columns = {f"_atom_site.{field}": data[f"_atom_site.{field}"] for field in fields}
    if len({len(values) for values in columns.values()}) != 1:
        raise ValueError(f"{path.name} has inconsistent atom columns")
    best: dict[tuple[str, int, str], tuple[float, str, str, int, str, np.ndarray]] = {}
    ligand: dict[str, np.ndarray] = {}
    for i in range(len(columns["_atom_site.group_PDB"])):
        if columns["_atom_site.pdbx_PDB_model_num"][i] != "1":
            continue
        group = columns["_atom_site.group_PDB"][i]
        atom = columns["_atom_site.label_atom_id"][i]
        comp = columns["_atom_site.label_comp_id"][i]
        chain = columns["_atom_site.label_asym_id"][i]
        if group == "HETATM" and comp == "ADP" and columns["_atom_site.auth_asym_id"][i] == "A" and atom in ADENINE_ATOMS:
            ligand[atom] = _coordinate(columns, i)
            continue
        if group != "ATOM" or chain not in chains or columns["_atom_site.type_symbol"][i] in ("H", "D"):
            continue
        try:
            auth_seq = int(columns["_atom_site.auth_seq_id"][i])
            label_seq = int(columns["_atom_site.label_seq_id"][i])
        except ValueError:
            continue
        occupancy = float(columns["_atom_site.occupancy"][i])
        alt = columns["_atom_site.label_alt_id"][i]
        key = (chain, auth_seq, atom)
        candidate = (occupancy, alt, comp, label_seq, columns["_atom_site.auth_asym_id"][i], _coordinate(columns, i))
        old = best.get(key)
        if old is None or (occupancy, alt == ".", alt == "A") > (old[0], old[1] == ".", old[1] == "A"):
            best[key] = candidate
    grouped: dict[tuple[str, int], list[tuple[str, tuple[float, str, str, int, str, np.ndarray]]]] = {}
    for (chain, pos, atom), record in best.items():
        grouped.setdefault((chain, pos), []).append((atom, record))
    sites: dict[str, dict[int, Site]] = {chain: {} for chain in chains}
    for (chain, pos), atom_records in grouped.items():
        identities = {(record[2], record[3], record[4]) for _, record in atom_records}
        if len(identities) != 1:
            raise ValueError(f"ambiguous identity/mapping at {path.name} {chain}:{pos}")
        comp, label_seq, auth_chain = identities.pop()
        letter = residue_letter(comp)
        if chain == "A" and (pos < 1 or pos > len(sequence) or sequence[pos - 1] != letter):
            raise ValueError(f"MYH7/P12883 identity mismatch at {path.name} {chain}:{pos}")
        atom_dict = {atom: record[5] for atom, record in atom_records}
        sites[chain][pos] = Site(
            residue=letter,
            auth_chain=auth_chain,
            auth_seq=pos,
            label_seq=label_seq,
            ca=atom_dict.get("CA"),
            atoms=np.asarray(list(atom_dict.values()), dtype=float),
        )
    return CifStructure(title=data.get("_struct.title", [""])[0], sites=sites, ligand_atoms=ligand)


def load_8act(path: Path, sequence: str) -> dict[str, dict[int, Site]]:
    model = PDBParser(QUIET=True).get_structure("8ACT", str(path))[0]
    required = {"A", "B", "C", "D"}
    if not required.issubset({chain.id for chain in model}):
        raise ValueError("8ACT is missing heavy or essential-light chains")
    sites: dict[str, dict[int, Site]] = {chain: {} for chain in required}
    for chain in required:
        for residue in model[chain]:
            if residue.id[0] != " ":
                continue
            pos = residue.id[1]
            letter = residue_letter(residue.resname)
            if chain in ("A", "B") and (pos < 1 or pos > len(sequence) or sequence[pos - 1] != letter):
                raise ValueError(f"8ACT/P12883 identity mismatch at {chain}:{pos}")
            atoms = np.asarray([atom.coord for atom in residue if atom.element not in ("H", "D")], dtype=float)
            ca = residue["CA"].coord if "CA" in residue else None
            sites[chain][pos] = Site(letter, chain, pos, pos, ca, atoms)
    return sites


def minimum_distance(a: np.ndarray, b: np.ndarray) -> float:
    if not len(a) or not len(b):
        raise ValueError("distance requires observed atoms")
    return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=2).min())


def align_adp_to_rigor(adp: dict[int, Site], rigor: dict[int, Site]) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    positions = [pos for start, end in ALIGNMENT_CORE for pos in range(start, end + 1)
                 if pos in adp and pos in rigor and adp[pos].ca is not None and rigor[pos].ca is not None]
    if len(positions) < 300:
        raise ValueError("insufficient shared motor-core C-alpha coverage for alignment")
    moving = np.asarray([adp[pos].ca for pos in positions], dtype=float)
    target = np.asarray([rigor[pos].ca for pos in positions], dtype=float)
    moving_center = moving.mean(axis=0)
    target_center = target.mean(axis=0)
    u, _, vt = np.linalg.svd((moving - moving_center).T @ (target - target_center))
    parity = np.linalg.det(u @ vt)
    rotation = u @ np.diag([1.0, 1.0, 1.0 if parity >= 0 else -1.0]) @ vt
    return moving_center, target_center, rotation


def _site_mapping(variant: str, position: int, wt: str, entry: str, chain: str, site: Site | None) -> dict[str, str | int]:
    return {
        "variant": variant, "canonical_position": position, "canonical_wt": wt,
        "structure": entry, "label_chain": chain,
        "author_chain": site.auth_chain if site else "",
        "label_seq_id": site.label_seq if site else "",
        "author_seq_id": site.auth_seq if site else "",
        "observed_wt": site.residue if site else "",
        "has_ca": "1" if site is not None and site.ca is not None else "0",
        "heavy_atom_count": len(site.atoms) if site else 0,
        "status": "observed" if site is not None and site.ca is not None else "missing_ca_or_residue",
    }


def _write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def generate_features(*, variants_path: Path, fasta_path: Path, assembly_path: Path,
                      adp_path: Path, rigor_path: Path, actin_path: Path, output_dir: Path) -> None:
    paths = {"fasta": fasta_path, "8ACT": assembly_path, "8EFE": adp_path, "8EFD": rigor_path, "8EFI": actin_path}
    hashes = {name: verify_input(name, path) for name, path in paths.items()}
    sequence = load_fasta(fasta_path)
    variants = load_variants(variants_path, sequence)
    adp = load_cif(adp_path, sequence, chains={"A"})
    rigor = load_cif(rigor_path, sequence, chains={"A"})
    actomyosin = load_cif(actin_path, sequence, chains={"A", "B", "C", "D", "E", "F"})
    folded = load_8act(assembly_path, sequence)
    if len(adp.ligand_atoms) < 8:
        raise ValueError("8EFE lacks sufficient deposited ADP adenine atoms")
    for chain in ("A", "B"):
        for landmark in RELAY_LANDMARKS:
            if landmark not in folded[chain] or folded[chain][landmark].ca is None:
                raise ValueError(f"8ACT missing relay landmark {chain}:{landmark}")
    center_adp, center_rigor, rotation = align_adp_to_rigor(adp.sites["A"], rigor.sites["A"])
    actin_atoms = np.concatenate([site.atoms for chain in "BCDEF" for site in actomyosin.sites[chain].values()])
    elc_atoms = {
        "A": np.concatenate([site.atoms for site in folded["D"].values()]),
        "B": np.concatenate([site.atoms for site in folded["C"].values()]),
    }
    feature_rows: list[dict] = []
    coverage_rows: list[dict] = []
    mapping_rows: list[dict] = []
    for variant, pos, wt in variants:
        row: dict[str, str | float] = {"variant": variant}
        selected = {
            ("8EFE", "A"): adp.sites["A"].get(pos),
            ("8EFD", "A"): rigor.sites["A"].get(pos),
            ("8EFI", "A"): actomyosin.sites["A"].get(pos),
            ("8ACT", "A"): folded["A"].get(pos),
            ("8ACT", "B"): folded["B"].get(pos),
        }
        for (entry, chain), site in selected.items():
            mapping_rows.append(_site_mapping(variant, pos, wt, entry, chain, site))
            if site is not None and site.residue != wt:
                raise ValueError(f"{variant} WT identity mismatch in {entry} chain {chain}")

        def record(feature: str, value: float | None, reason: str = "") -> None:
            row[feature] = "" if value is None else round(value, 9)
            coverage_rows.append({"variant": variant, "feature": feature,
                                  "available": "1" if value is not None else "0", "missing_reason": reason})

        adp_site = selected[("8EFE", "A")]
        rigor_site = selected[("8EFD", "A")]
        if adp_site is None or rigor_site is None or adp_site.ca is None or rigor_site.ca is None:
            record(FEATURES[0], None, "site_not_observed_in_8EFE_or_8EFD")
        else:
            aligned = (adp_site.ca - center_adp) @ rotation + center_rigor
            record(FEATURES[0], float(np.linalg.norm(aligned - rigor_site.ca)))
        if adp_site is None or adp_site.ca is None:
            record(FEATURES[1], None, "site_not_observed_in_8EFE")
        else:
            record(FEATURES[1], minimum_distance(adp_site.ca[None, :], np.asarray(list(adp.ligand_atoms.values()))))
        actin_site = selected[("8EFI", "A")]
        if actin_site is None or not len(actin_site.atoms):
            record(FEATURES[2], None, "site_not_observed_in_8EFI")
        else:
            record(FEATURES[2], minimum_distance(actin_site.atoms, actin_atoms))
        a, b = selected[("8ACT", "A")], selected[("8ACT", "B")]
        if a is None or b is None or a.ca is None or b.ca is None:
            record(FEATURES[3], None, "site_not_observed_on_both_8ACT_heavy_chains")
            record(FEATURES[4], None, "site_not_observed_on_both_8ACT_heavy_chains")
        else:
            distances = []
            elc_distances = []
            for chain, site in (("A", a), ("B", b)):
                landmarks = np.asarray([folded[chain][p].ca for p in RELAY_LANDMARKS])
                distances.append(minimum_distance(site.ca[None, :], landmarks))
                elc_distances.append(minimum_distance(site.atoms, elc_atoms[chain]))
            record(FEATURES[3], float(np.mean(distances)))
            record(FEATURES[4], float(np.mean(elc_distances)))
        feature_rows.append(row)
    output_dir.mkdir(parents=True, exist_ok=True)
    feature_path = output_dir / "features-v1.csv"
    coverage_path = output_dir / "feature-coverage-v1.csv"
    mapping_path = output_dir / "site-mapping-v1.csv"
    _write_csv(feature_path, feature_rows, ["variant", *FEATURES])
    _write_csv(coverage_path, coverage_rows, ["variant", "feature", "available", "missing_reason"])
    _write_csv(mapping_path, mapping_rows, ["variant", "canonical_position", "canonical_wt", "structure", "label_chain", "author_chain", "label_seq_id", "author_seq_id", "observed_wt", "has_ca", "heavy_atom_count", "status"])
    manifest = {
        "schema_version": 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "feature_generation_uses_outcomes": False,
        "variant_registry_sha256": sha256(variants_path),
        "sequence_accession": "P12883",
        "normalized_sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "features": {feature: {"unit": "angstrom", "definition": "docs/myh7-correlation-feature-literature.md"} for feature in FEATURES},
        "structures": {
            entry: {
                "url": STRUCTURE_URLS[entry],
                "sha256": hashes[entry],
                "sha256_applies_to": "decompressed_PDB_input" if entry == "8ACT" else "downloaded_mmCIF_input",
                "download_sha256": "84469940f04dd8428a98e76c553ab519c1d45aef2beda176b8813d3b120e797b" if entry == "8ACT" else hashes[entry],
                "transform": "gzip_decompress" if entry == "8ACT" else "none",
                "input_kind": "deposited_WT_coordinates",
            }
            for entry in STRUCTURE_URLS
        },
        "fasta_sha256": hashes["fasta"],
        "mapping": {"8EFE": "label A/auth A human MYH7", "8EFD": "label A/auth A human MYH7", "8EFI": "label A/auth M human MYH7; actin label B-F", "8ACT": "auth A/B human MYH7; spatially paired ELC A-D/B-C"},
        "adenine_atoms_expected": sorted(ADENINE_ATOMS),
        "adenine_atoms_observed": sorted(adp.ligand_atoms),
        "alignment_core": [list(bounds) for bounds in ALIGNMENT_CORE],
        "relay_landmarks": list(RELAY_LANDMARKS),
        "output_sha256": {path.name: sha256(path) for path in (feature_path, coverage_path, mapping_path)},
        "tool_versions": {"biopython": biopython_version, "numpy": np.__version__},
    }
    (output_dir / "feature-manifest-v1.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
