"""Failing-first contract tests for 8ACT structure features (#28 Task 3)."""

from __future__ import annotations

import gzip
import hashlib
import io
import os
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest
from Bio.PDB import Atom, Chain, Model, Residue, Structure

from lsar_na_rel.structure import (
    ASSEMBLY_URL,
    CONTACT_CUTOFF_ANGSTROM,
    DECOMPRESSED_PDB_SHA256,
    DEPOSITED_HELIX_ONSETS,
    DEPOSITED_HELIX_RANGES,
    GZIP_SHA256,
    MOTOR_LEVER_RANGES,
    PROXIMAL_S2_RANGES,
    StructureResult,
    _download_url_to_path,
    ensure_8act_assembly,
    generate_structure_features,
    heavy_atom_contact,
    ihm_partner_contacts,
    load_8act_model,
    normalize_residue_identity,
    residues_in_heavy_contact,
    validate_wt_sites_ab,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
WORKFLOW_ROOT = REPO_ROOT / "workflows" / "lsar-na-rel"
CACHE_DIR = WORKFLOW_ROOT / "cache"

# Frozen primary Talk-v0 labels (order matches canonical S1 primary table).
# Kept as a local tuple so structure tests never load Na_rel / ln(Na_rel).
PRIMARY_VARIANTS: tuple[str, ...] = (
    "Y115H",
    "E497D",
    "R249Q",
    "H251N",
    "D382Y",
    "I457T",
    "R719W",
    "P710R",
    "D778V",
    "L781P",
    "S782N",
    "A797T",
    "F834L",
    "R403Q",
    "R663H",
    "G768R",
)

# Pinned RSA (A, B, mean) under mkdssp 4.6.1 / Biopython 1.88 / Sander.
PINNED_RSA = {
    115: (0.027027027, 0.009009009, 0.018018018),
    497: (0.144329897, 0.030927835, 0.087628866),
    249: (0.153225806, 0.258064516, 0.205645161),
    251: (0.114130435, 0.076086957, 0.095108696),
    382: (0.073619632, 0.392638037, 0.233128834),
    457: (0.0, 0.0, 0.0),
    719: (0.165322581, 0.157258065, 0.161290323),
    710: (0.139705882, 0.080882353, 0.110294118),
    778: (0.153374233, 0.478527607, 0.315950920),
    781: (0.0, 0.304878049, 0.152439024),
    782: (0.223076923, 0.461538462, 0.342307692),
    797: (0.037735849, 0.018867925, 0.028301887),
    834: (0.203045685, 0.218274112, 0.210659898),
    403: (0.137096774, 0.229838710, 0.183467742),
    663: (0.588709677, 0.629032258, 0.608870968),
    768: (0.5, 0.726190476, 0.613095238),
}

BH_FH_POSITIONS = frozenset({382, 403})
LIGHT_CHAIN_ONLY_POSITIONS = frozenset({778, 781, 782, 797, 834})


@pytest.fixture(scope="module")
def primary_variants() -> tuple[str, ...]:
    return PRIMARY_VARIANTS


@pytest.fixture(scope="module")
def assembly_manifest():
    return ensure_8act_assembly(CACHE_DIR)


@pytest.fixture(scope="module")
def assembly_model(assembly_manifest):
    return load_8act_model(assembly_manifest.pdb_path)


@pytest.fixture(scope="module")
def structure_result(primary_variants, assembly_manifest) -> StructureResult:
    return generate_structure_features(
        primary_variants,
        assembly_pdb_path=assembly_manifest.pdb_path,
        mkdssp="mkdssp",
    )


def _add_atom(residue: Residue, name: str, coord, element: str) -> Atom:
    atom = Atom.Atom(
        name,
        np.asarray(coord, dtype=float),
        20.0,
        1.0,
        " ",
        name,
        len(list(residue.get_atoms())) + 1,
        element,
    )
    residue.add(atom)
    return atom


def _synthetic_residue(chain_id: str, resseq: int, resname: str, ca_coord) -> Residue:
    residue = Residue.Residue((" ", resseq, " "), resname, " ")
    _add_atom(residue, "CA", ca_coord, "C")
    return residue


def _build_synthetic_model(residues: list[tuple[str, int, str, tuple[float, float, float]]]) -> Model:
    structure = Structure.Structure("synth")
    model = Model.Model(0)
    structure.add(model)
    chains: dict[str, Chain] = {}
    for chain_id, resseq, resname, coord in residues:
        if chain_id not in chains:
            chain = Chain.Chain(chain_id)
            chains[chain_id] = chain
            model.add(chain)
        chains[chain_id].add(_synthetic_residue(chain_id, resseq, resname, coord))
    return model


def _sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def test_checksummed_assembly_manifest_matches_pinned_hashes(assembly_manifest):
    assert assembly_manifest.url == ASSEMBLY_URL
    assert assembly_manifest.gzip_sha256 == GZIP_SHA256
    assert assembly_manifest.pdb_sha256 == DECOMPRESSED_PDB_SHA256
    assert assembly_manifest.gzip_sha256 == (
        "84469940f04dd8428a98e76c553ab519c1d45aef2beda176b8813d3b120e797b"
    )
    assert assembly_manifest.pdb_sha256 == (
        "ef548f824d10986e313fa41910fa10d5d3925bc04ae3350f58214a2639d3d6a9"
    )
    assert assembly_manifest.pdb_path.is_file()
    assert assembly_manifest.gzip_path.is_file()


def test_chains_a_through_f_parse(assembly_model):
    chain_ids = set(chain.id for chain in assembly_model)
    assert chain_ids == {"A", "B", "C", "D", "E", "F"}


def test_load_8act_model_accepts_any_order_of_required_chains(tmp_path):
    # Minimal PDB with chains in non-canonical order F..A.
    lines = []
    serial = 1
    for chain_id in ("F", "E", "D", "C", "B", "A"):
        x = float(serial)
        # Exact PDB fixed-column ATOM record.
        lines.append(
            f"ATOM  {serial:5d}  CA  ALA {chain_id}   1    "
            f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00 20.00           C"
        )
        serial += 1
    lines.append("END")
    pdb_path = tmp_path / "reordered.pdb"
    pdb_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    model = load_8act_model(pdb_path)
    assert {chain.id for chain in model} == {"A", "B", "C", "D", "E", "F"}


def test_all_16_positions_match_expected_wt_in_a_and_b(assembly_model, primary_variants):
    sites = validate_wt_sites_ab(assembly_model, primary_variants)
    assert len(sites) == 16
    for label, site in sites.items():
        assert site["A"] == site["claimed_wt"]
        assert site["B"] == site["claimed_wt"]
        assert site["covered_a"] is True
        assert site["covered_b"] is True
    # Position-independent M3L→K normalization for identity.
    assert normalize_residue_identity("M3L") == "K"
    assert normalize_residue_identity("LYS") == "K"


def test_m3l_at_806_accepts_claimed_lysine_identity_on_both_chains():
    """Explicit 806 coverage: M3L maps to K for identity only.

    May pass immediately under the existing position-independent M3L→K rule.
    """
    model = _build_synthetic_model(
        [
            ("A", 806, "M3L", (0.0, 0.0, 0.0)),
            ("B", 806, "M3L", (10.0, 0.0, 0.0)),
        ]
    )
    sites = validate_wt_sites_ab(model, ["K806A"])
    assert sites["K806A"]["A"] == "K"
    assert sites["K806A"]["B"] == "K"
    assert sites["K806A"]["claimed_wt"] == "K"


def test_dssp_invoked_with_sander_acc_array(assembly_manifest, primary_variants):
    sentinel = object()

    def _fake_dssp(model, in_file, dssp="dssp", acc_array="Sander", file_type=""):
        assert Path(in_file) == assembly_manifest.pdb_path
        assert dssp == "mkdssp"
        assert acc_array == "Sander"
        return sentinel

    with patch("lsar_na_rel.structure.DSSP", side_effect=_fake_dssp) as mock_dssp:
        with patch(
            "lsar_na_rel.structure._rsa_from_dssp",
            return_value={
                ("A", 115): 0.1,
                ("B", 115): 0.2,
            },
        ):
            generate_structure_features(
                ["Y115H"],
                assembly_pdb_path=assembly_manifest.pdb_path,
                mkdssp="mkdssp",
            )
    assert mock_dssp.called


def test_generate_rejects_corrupt_pdb_before_dssp(tmp_path, assembly_manifest):
    corrupt = tmp_path / "corrupt-8ACT.pdb"
    corrupt.write_bytes(b"HEADER CORRUPT SUBSTITUTED PDB\nEND\n")
    assert _sha256_bytes(corrupt.read_bytes()) != DECOMPRESSED_PDB_SHA256

    with patch("lsar_na_rel.structure.DSSP") as mock_dssp:
        with pytest.raises(ValueError, match="SHA-256|checksum|mismatch"):
            generate_structure_features(
                ["Y115H"],
                assembly_pdb_path=corrupt,
                mkdssp="mkdssp",
            )
    assert mock_dssp.call_count == 0


def test_generate_provenance_records_observed_verified_pdb_hash(structure_result):
    provenance = structure_result.provenance
    assert provenance["pdb_sha256"] == DECOMPRESSED_PDB_SHA256
    assert provenance["pdb_sha256_verified"] is True
    # Do not claim a gzip hash was verified on the PDB-only generation path.
    assert "gzip_sha256" not in provenance or provenance.get("gzip_sha256_verified") is not True
    assert provenance.get("pinned_pdb_sha256") == DECOMPRESSED_PDB_SHA256
    assert provenance.get("pinned_gzip_sha256") == GZIP_SHA256
    assert provenance.get("assembly_url") == ASSEMBLY_URL


def test_pinned_rsa_mean_ab_matches_sander_normalization(structure_result, primary_variants):
    features = structure_result.features.set_index("variant")
    sites = structure_result.site_provenance.set_index("variant")
    for label in primary_variants:
        position = int("".join(ch for ch in label if ch.isdigit()))
        rsa_a, rsa_b, mean = PINNED_RSA[position]
        assert sites.loc[label, "wt_site_rsa_8act_a"] == pytest.approx(rsa_a, abs=1e-9)
        assert sites.loc[label, "wt_site_rsa_8act_b"] == pytest.approx(rsa_b, abs=1e-9)
        assert features.loc[label, "wt_site_rsa_8act_mean_ab"] == pytest.approx(mean, abs=1e-9)
        assert sites.loc[label, "wt_site_rsa_8act_bh"] == pytest.approx(rsa_a, abs=1e-9)
        assert sites.loc[label, "wt_site_rsa_8act_fh"] == pytest.approx(rsa_b, abs=1e-9)
        assert sites.loc[label, "coverage_a"] is True or sites.loc[label, "coverage_a"] == 1
        assert sites.loc[label, "coverage_b"] is True or sites.loc[label, "coverage_b"] == 1


def test_heavy_atom_contact_boundary_ignores_hydrogens():
    # Exactly 4.5 Å between heavy atoms → True; greater → False.
    assert heavy_atom_contact((0.0, 0.0, 0.0), (4.5, 0.0, 0.0), cutoff=4.5) is True
    assert heavy_atom_contact((0.0, 0.0, 0.0), (4.5000001, 0.0, 0.0), cutoff=4.5) is False

    residue_a = Residue.Residue((" ", 1, " "), "ALA", " ")
    residue_b = Residue.Residue((" ", 2, " "), "ALA", " ")
    _add_atom(residue_a, "CA", (0.0, 0.0, 0.0), "C")
    _add_atom(residue_a, "H", (0.0, 0.0, 0.1), "H")
    # Hydrogen alone would be within 1 Å, but must be ignored; heavy atoms at 5 Å.
    _add_atom(residue_b, "CA", (5.0, 0.0, 0.0), "C")
    _add_atom(residue_b, "H", (0.5, 0.0, 0.0), "H")
    assert residues_in_heavy_contact(residue_a, residue_b, cutoff=4.5) is False


def test_blank_element_digit_prefixed_hydrogen_ignored_in_contacts():
    import warnings

    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    residue_a = Residue.Residue((" ", 1, " "), "ALA", " ")
    residue_b = Residue.Residue((" ", 2, " "), "ALA", " ")
    _add_atom(residue_a, "CA", (0.0, 0.0, 0.0), "C")
    # Force blank element after construction: Biopython may auto-fill from the
    # atom name, but deposited PDBs can leave digit-prefixed hydrogens blank.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        h_a = _add_atom(residue_a, "1HB", (0.2, 0.0, 0.0), "")
        h_b = _add_atom(residue_b, "1HB", (0.5, 0.0, 0.0), "")
    h_a.element = " "
    _add_atom(residue_b, "CA", (5.0, 0.0, 0.0), "C")
    h_b.element = " "
    assert residues_in_heavy_contact(residue_a, residue_b, cutoff=4.5) is False


def test_deuterium_element_and_digit_prefixed_blank_name_ignored_in_contacts():
    import warnings

    from Bio.PDB.PDBExceptions import PDBConstructionWarning

    # Explicit element D near a far heavy atom must not create a contact.
    residue_a = Residue.Residue((" ", 1, " "), "ALA", " ")
    residue_b = Residue.Residue((" ", 2, " "), "ALA", " ")
    _add_atom(residue_a, "CA", (0.0, 0.0, 0.0), "C")
    _add_atom(residue_a, "D", (0.2, 0.0, 0.0), "D")
    _add_atom(residue_b, "CA", (5.0, 0.0, 0.0), "C")
    _add_atom(residue_b, "D", (0.4, 0.0, 0.0), "D")
    assert residues_in_heavy_contact(residue_a, residue_b, cutoff=4.5) is False

    # Digit-prefixed deuterium name with blank element (e.g. 1D).
    residue_c = Residue.Residue((" ", 3, " "), "ALA", " ")
    residue_d = Residue.Residue((" ", 4, " "), "ALA", " ")
    _add_atom(residue_c, "CA", (0.0, 0.0, 0.0), "C")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=PDBConstructionWarning)
        d_c = _add_atom(residue_c, "1D", (0.2, 0.0, 0.0), "")
        d_d = _add_atom(residue_d, "1D", (0.5, 0.0, 0.0), "")
    d_c.element = " "
    _add_atom(residue_d, "CA", (5.0, 0.0, 0.0), "C")
    d_d.element = " "
    assert residues_in_heavy_contact(residue_c, residue_d, cutoff=4.5) is False


def test_synthetic_ihm_partner_classification_bh_fh_bh_s2_and_light_only():
    # BH (A motor) contacts FH (B motor) at 4.0 Å → BH–FH.
    model_bh_fh = _build_synthetic_model(
        [
            ("A", 100, "ALA", (0.0, 0.0, 0.0)),
            ("B", 200, "ALA", (4.0, 0.0, 0.0)),
        ]
    )
    flags = ihm_partner_contacts(model_bh_fh, position=100)
    assert flags["ihm_contact_bh_fh"] == 1
    assert flags["ihm_contact_bh_s2"] == 0
    assert flags["ihm_interface_flag_8act"] == 1

    # BH (A motor) contacts partner proximal S2 (B >=838) → BH–S2.
    model_bh_s2 = _build_synthetic_model(
        [
            ("A", 100, "ALA", (0.0, 0.0, 0.0)),
            ("B", 850, "ALA", (4.0, 0.0, 0.0)),
        ]
    )
    flags_s2 = ihm_partner_contacts(model_bh_s2, position=100)
    assert flags_s2["ihm_contact_bh_fh"] == 0
    assert flags_s2["ihm_contact_bh_s2"] == 1
    assert flags_s2["ihm_interface_flag_8act"] == 1

    # Light-chain-only contact (C) must not set either subtype or main flag.
    model_light = _build_synthetic_model(
        [
            ("A", 100, "ALA", (0.0, 0.0, 0.0)),
            ("C", 50, "ALA", (3.0, 0.0, 0.0)),
        ]
    )
    flags_light = ihm_partner_contacts(model_light, position=100)
    assert flags_light["ihm_contact_bh_fh"] == 0
    assert flags_light["ihm_contact_bh_s2"] == 0
    assert flags_light["ihm_interface_flag_8act"] == 0


def test_s1_s2_sequence_boundary_837_motor_vs_838_proximal_s2():
    """Exact S1/S2 cut: partner 837 is BH–FH; partner 838 is BH–S2."""
    model_motor = _build_synthetic_model(
        [
            ("A", 100, "ALA", (0.0, 0.0, 0.0)),
            ("B", 837, "ALA", (4.0, 0.0, 0.0)),
        ]
    )
    flags_motor = ihm_partner_contacts(model_motor, position=100)
    assert flags_motor["ihm_contact_bh_fh"] == 1
    assert flags_motor["ihm_contact_bh_s2"] == 0
    assert flags_motor["ihm_interface_flag_8act"] == 1

    model_s2 = _build_synthetic_model(
        [
            ("A", 100, "ALA", (0.0, 0.0, 0.0)),
            ("B", 838, "ALA", (4.0, 0.0, 0.0)),
        ]
    )
    flags_s2 = ihm_partner_contacts(model_s2, position=100)
    assert flags_s2["ihm_contact_bh_fh"] == 0
    assert flags_s2["ihm_contact_bh_s2"] == 1
    assert flags_s2["ihm_interface_flag_8act"] == 1


def test_subtype_bits_and_main_or_on_real_assembly(structure_result, primary_variants):
    sites = structure_result.site_provenance.set_index("variant")
    features = structure_result.features.set_index("variant")
    for label in primary_variants:
        position = int("".join(ch for ch in label if ch.isdigit()))
        bh_fh = int(sites.loc[label, "ihm_contact_bh_fh"])
        bh_s2 = int(sites.loc[label, "ihm_contact_bh_s2"])
        main = int(features.loc[label, "ihm_interface_flag_8act"])
        assert main == (1 if (bh_fh or bh_s2) else 0)
        assert bh_s2 == 0
        if position in BH_FH_POSITIONS:
            assert bh_fh == 1
            assert main == 1
        else:
            assert bh_fh == 0
            assert main == 0
        if position in LIGHT_CHAIN_ONLY_POSITIONS:
            assert main == 0


def test_features_contain_only_rsa_mean_and_ihm_flag(structure_result):
    assert list(structure_result.features.columns) == [
        "variant",
        "wt_site_rsa_8act_mean_ab",
        "ihm_interface_flag_8act",
    ]
    forbidden = {
        "delta_volume",
        "volume",
        "delta_hydropathy",
        "hydropathy",
        "contact_count",
        "contact_count_8act_mean_ab",
        "foldx_ddg",
        "foldx_ddg_bh_fh",
        "foldx_ddg_mean_ab",
        "secondary_structure",
        "domain",
        "wt_site_rsa_8act_a",
        "wt_site_rsa_8act_b",
        "ihm_contact_bh_fh",
        "ihm_contact_bh_s2",
    }
    assert set(structure_result.features.columns).isdisjoint(forbidden)


def test_structure_does_not_load_target_or_labels(structure_result, primary_variants, assembly_manifest):
    assert isinstance(structure_result, StructureResult)
    assert "y" not in structure_result.features.columns
    assert "Na_rel" not in structure_result.features.columns
    assert "ln(Na_rel)" not in structure_result.features.columns
    assert "source_id" not in structure_result.features.columns
    assert not hasattr(structure_result, "target")

    def _forbid_read_csv(*_args, **_kwargs):
        raise AssertionError("structure generation must not call pd.read_csv")

    with patch("pandas.read_csv", side_effect=_forbid_read_csv):
        with patch("lsar_na_rel.structure.pd.read_csv", side_effect=_forbid_read_csv):
            generate_structure_features(
                primary_variants,
                assembly_pdb_path=assembly_manifest.pdb_path,
                mkdssp="mkdssp",
            )


def test_ensure_8act_assembly_downloads_from_file_url_when_gzip_absent(tmp_path):
    pdb_bytes = b"HEADER TEST ASSEMBLY\nATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n"
    gzip_buf = io.BytesIO()
    with gzip.GzipFile(fileobj=gzip_buf, mode="wb") as handle:
        handle.write(pdb_bytes)
    gzip_bytes = gzip_buf.getvalue()
    gzip_hash = _sha256_bytes(gzip_bytes)
    pdb_hash = _sha256_bytes(pdb_bytes)

    source = tmp_path / "source.pdb1.gz"
    source.write_bytes(gzip_bytes)
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    url = source.resolve().as_uri()

    manifest = ensure_8act_assembly(
        cache_dir,
        url=url,
        gzip_sha256=gzip_hash,
        pdb_sha256=pdb_hash,
        gzip_name="toy.pdb1.gz",
        pdb_name="toy-assembly1.pdb",
    )
    assert manifest.gzip_path.is_file()
    assert manifest.pdb_path.is_file()
    assert manifest.gzip_sha256 == gzip_hash
    assert manifest.pdb_sha256 == pdb_hash
    assert manifest.pdb_path.read_bytes() == pdb_bytes


def test_ensure_8act_assembly_rejects_bad_download_and_cleans_partial(tmp_path):
    pdb_bytes = b"HEADER GOOD\nEND\n"
    gzip_buf = io.BytesIO()
    with gzip.GzipFile(fileobj=gzip_buf, mode="wb") as handle:
        handle.write(pdb_bytes)
    gzip_bytes = gzip_buf.getvalue()
    source = tmp_path / "source.pdb1.gz"
    source.write_bytes(gzip_bytes)
    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    url = source.resolve().as_uri()

    with pytest.raises(ValueError, match="SHA-256|checksum|mismatch"):
        ensure_8act_assembly(
            cache_dir,
            url=url,
            gzip_sha256="0" * 64,
            pdb_sha256=_sha256_bytes(pdb_bytes),
            gzip_name="bad.pdb1.gz",
            pdb_name="bad-assembly1.pdb",
        )
    assert not (cache_dir / "bad.pdb1.gz").exists()
    assert not (cache_dir / "bad-assembly1.pdb").exists()
    # No leftover temp partials.
    assert list(cache_dir.glob("*.partial")) == []
    assert list(cache_dir.glob("*.tmp")) == []


def test_failed_download_removes_unverified_final_target(tmp_path, monkeypatch):
    pdb_bytes = b"HEADER DOWNLOAD FAIL PATH\nEND\n"
    gzip_buf = io.BytesIO()
    with gzip.GzipFile(fileobj=gzip_buf, mode="wb") as handle:
        handle.write(pdb_bytes)
    source = tmp_path / "source.pdb1.gz"
    source.write_bytes(gzip_buf.getvalue())
    target = tmp_path / "cache" / "out.pdb1.gz"
    target.parent.mkdir()
    url = source.resolve().as_uri()

    real_replace = os.replace

    def _replace_then_fail(src, dst):
        real_replace(src, dst)
        raise OSError("simulated post-replace failure")

    monkeypatch.setattr(os, "replace", _replace_then_fail)
    with pytest.raises(OSError, match="simulated post-replace failure"):
        _download_url_to_path(url, target)
    assert not target.exists()
    assert list(target.parent.glob("*.partial")) == []


def test_ensure_8act_assembly_repairs_corrupt_pdb_from_verified_gzip(tmp_path):
    pdb_bytes = b"HEADER GOOD ASSEMBLY\nATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 20.00           C  \nEND\n"
    gzip_buf = io.BytesIO()
    with gzip.GzipFile(fileobj=gzip_buf, mode="wb") as handle:
        handle.write(pdb_bytes)
    gzip_bytes = gzip_buf.getvalue()
    gzip_hash = _sha256_bytes(gzip_bytes)
    pdb_hash = _sha256_bytes(pdb_bytes)

    cache_dir = tmp_path / "cache"
    cache_dir.mkdir()
    gzip_path = cache_dir / "toy.pdb1.gz"
    pdb_path = cache_dir / "toy-assembly1.pdb"
    gzip_path.write_bytes(gzip_bytes)
    pdb_path.write_bytes(b"HEADER CORRUPT LOCAL PDB\nEND\n")
    assert _sha256_bytes(pdb_path.read_bytes()) != pdb_hash

    manifest = ensure_8act_assembly(
        cache_dir,
        url="file:///unused-because-gzip-present",
        gzip_sha256=gzip_hash,
        pdb_sha256=pdb_hash,
        gzip_name="toy.pdb1.gz",
        pdb_name="toy-assembly1.pdb",
    )
    assert manifest.pdb_sha256 == pdb_hash
    assert pdb_path.read_bytes() == pdb_bytes


def test_download_urlopen_uses_finite_named_timeout(tmp_path, monkeypatch):
    from lsar_na_rel import structure as structure_mod

    assert hasattr(structure_mod, "DOWNLOAD_TIMEOUT_SECONDS")
    timeout = structure_mod.DOWNLOAD_TIMEOUT_SECONDS
    assert isinstance(timeout, (int, float))
    assert timeout > 0

    pdb_bytes = b"HEADER TIMEOUT\nEND\n"
    gzip_buf = io.BytesIO()
    with gzip.GzipFile(fileobj=gzip_buf, mode="wb") as handle:
        handle.write(pdb_bytes)
    source = tmp_path / "source.pdb1.gz"
    source.write_bytes(gzip_buf.getvalue())
    target = tmp_path / "out.pdb1.gz"
    url = source.resolve().as_uri()

    recorded: dict[str, object] = {}
    real_urlopen = __import__("urllib.request").request.urlopen

    def _tracking_urlopen(resource, *args, **kwargs):
        recorded["timeout"] = kwargs.get("timeout")
        return real_urlopen(resource, timeout=kwargs.get("timeout"))

    monkeypatch.setattr(structure_mod.urllib.request, "urlopen", _tracking_urlopen)
    _download_url_to_path(url, target)
    assert recorded["timeout"] == timeout
    assert target.is_file()


def test_frozen_partner_ranges_and_cutoff_are_exported():
    # Symmetric S1/S2 sequence cut (Pro838 = S1 C-terminus / S2 start), bounded
    # by modeled 8ACT residues. Deposited helix onsets A841/B846 are provenance.
    assert MOTOR_LEVER_RANGES == {"A": (3, 837), "B": (3, 837)}
    assert PROXIMAL_S2_RANGES == {"A": (838, 906), "B": (838, 898)}
    assert DEPOSITED_HELIX_ONSETS == {"A": 841, "B": 846}
    assert DEPOSITED_HELIX_RANGES == {"A": (841, 906), "B": (846, 898)}
    assert CONTACT_CUTOFF_ANGSTROM == 4.5


def test_heavy_atom_residue_contact_at_exact_boundary():
    residue_a = Residue.Residue((" ", 1, " "), "ALA", " ")
    residue_b = Residue.Residue((" ", 2, " "), "ALA", " ")
    _add_atom(residue_a, "CA", (0.0, 0.0, 0.0), "C")
    _add_atom(residue_b, "CA", (4.5, 0.0, 0.0), "C")

    assert residues_in_heavy_contact(residue_a, residue_b, cutoff=4.5) is True
    residue_far = Residue.Residue((" ", 3, " "), "ALA", " ")
    _add_atom(residue_far, "CA", (4.5000001, 0.0, 0.0), "C")
    assert residues_in_heavy_contact(residue_a, residue_far, cutoff=4.5) is False
    # Hydrogens within cutoff must not create a contact when heavy atoms are far.
    residue_h = Residue.Residue((" ", 4, " "), "ALA", " ")
    _add_atom(residue_h, "CA", (5.0, 0.0, 0.0), "C")
    _add_atom(residue_h, "H1", (0.1, 0.0, 0.0), "H")
    _add_atom(residue_a, "H2", (0.0, 0.1, 0.0), "H")
    assert residues_in_heavy_contact(residue_a, residue_h, cutoff=4.5) is False
