"""Golden tests for Phase 2 transforms.

Each test uses dict-literal fixtures and a frozen clock.
No database, no filesystem, no IO.
"""
import datetime
import copy

# Import the transforms and helpers from lw.py (same directory)
import sys, os
sys.path.insert(0, os.path.dirname(__file__))
from lw_refactoring import (
    transform_exp_complete,
    transform_exp_show,
    effects_ok,
    effects_fail,
    lookup_experiment,
)


# ============================================================
# FIXTURES
# ============================================================

FROZEN_CLOCK = {
    "today": datetime.date(2025, 6, 10),
    "now": datetime.datetime(2025, 6, 10, 14, 30, 0),
}

MINIMAL_STORE = {
    "experiments": [
        {
            "id": "LM-0001",
            "type": "purification",
            "title": "test-pur",
            "shortname": "test-pur",
            "date_started": "2025-06-01",
            "date_completed": None,
            "status": "active",
            "folder_name": "20250601_LM-0001_pur_test-pur",
            "protocol_id": None,
            "notes": None,
            "created_at": "2025-06-01T00:00:00",
        },
    ],
    "strains": [],
    "experiment_strains": [],
    "experiment_links": [],
    "files": [],
    "purification_details": [],
    "assay_details": [],
    "figure_panels": [],
    "samples": [],
    "stage_expectations": [],
}

COMPLETE_EXPERIMENT_STORE = {
    "experiments": [
        {
            "id": "LM-0002",
            "type": "gelshift",
            "title": "gel-test",
            "shortname": "gel-test",
            "date_started": "2025-05-20",
            "date_completed": "2025-06-05",
            "status": "complete",
            "folder_name": "20250520_LM-0002_gs_gel-test",
            "protocol_id": "P-001",
            "notes": "[2025-05-20] initial notes",
            "created_at": "2025-05-20T00:00:00",
        },
    ],
    "strains": [
        {
            "id": "LY456",
            "genotype": "MATa orc4-R267A",
            "label": "orc4-mut",
            "parent": None,
            "construction": None,
            "selection_marker": "KanMX",
            "verification": None,
            "storage_location": None,
            "notes": None,
        },
    ],
    "experiment_strains": [
        {"experiment_id": "LM-0002", "strain_id": "LY456", "role": "test"},
    ],
    "experiment_links": [
        {
            "from_experiment": "LM-0002",
            "to_experiment": "LM-0001",
            "relationship": "uses_prep_from",
        },
        {
            "from_experiment": "LM-0003",
            "to_experiment": "LM-0002",
            "relationship": "replicate_of",
        },
    ],
    "files": [
        {
            "id": 1,
            "experiment_id": "LM-0002",
            "filename": "experiments/20250520_LM-0002_gs_gel-test/20250521_LM-0002_gel-01.tiff",
            "file_type": "gel",
            "instrument": None,
            "description": "coomassie",
            "source_file_id": None,
            "created_at": "2025-05-21T00:00:00",
        },
    ],
    "purification_details": [],
    "assay_details": [
        {
            "experiment_id": "LM-0002",
            "protein_prep": "LM-0001",
            "dna_template": "pUC19",
            "dna_prep_method": None,
            "conditions": "37C 1hr",
        },
    ],
    "figure_panels": [],
    "samples": [
        {
            "id": 1,
            "experiment_id": "LM-0002",
            "position": 1,
            "sample_type": "factorial",
            "strain_id": "LY456",
            "dna": "pUC19",
            "label": "orc4-mut",
            "conditions": '{"strain": "LY456", "dna": "pUC19"}',
            "created_at": "2025-05-20T00:00:00",
        },
        {
            "id": 2,
            "experiment_id": "LM-0002",
            "position": 2,
            "sample_type": "extra",
            "strain_id": None,
            "dna": None,
            "label": "no protein",
            "conditions": "{}",
            "created_at": "2025-05-20T00:00:00",
        },
    ],
    "stage_expectations": [],
}


# ============================================================
# HELPERS
# ============================================================

def assert_ok(effects):
    assert effects["ok"], f"Expected ok=True, got stderr={effects['stderr']}"
    assert effects["exit_code"] == 0


def assert_fail(effects, exit_code=1):
    assert not effects["ok"], f"Expected ok=False, got stdout={effects['stdout']}"
    assert effects["exit_code"] == exit_code
    assert effects["store"] is None


def stdout_text(effects):
    return "\n".join(effects["stdout"])


# ============================================================
# TESTS: transform_exp_complete
# ============================================================

def test_complete_happy_path():
    """Active experiment -> complete. Store mutated, stdout matches."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_complete(store, {"exp_id": "LM-0001"}, FROZEN_CLOCK)
    assert_ok(effects)
    assert effects["store"] is not None
    # verify store mutation
    exp = lookup_experiment(effects["store"], "LM-0001")
    assert exp["status"] == "complete"
    assert exp["date_completed"] == "2025-06-10"
    # verify stdout
    text = stdout_text(effects)
    assert "Completed: LM-0001" in text
    assert "Date: 2025-06-10" in text
    assert "Duration: 9d (started 2025-06-01)" in text


def test_complete_already_complete():
    """Already-complete experiment -> no-op, exit_code 0, no store mutation."""
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_exp_complete(store, {"exp_id": "LM-0002"}, FROZEN_CLOCK)
    assert_ok(effects)
    assert effects["store"] is None  # no mutation
    text = stdout_text(effects)
    assert "already marked complete" in text


def test_complete_nonexistent():
    """Non-existent experiment -> error effects."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_complete(store, {"exp_id": "LM-9999"}, FROZEN_CLOCK)
    assert_fail(effects)
    assert "No experiment found" in effects["stderr"][0]


def test_complete_invalid_id():
    """Invalid ID -> error effects."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_complete(store, {"exp_id": "garbage"}, FROZEN_CLOCK)
    assert_fail(effects)
    assert "Invalid experiment ID" in effects["stderr"][0]


def test_complete_bare_number():
    """Bare number '1' normalizes to LM-0001."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_complete(store, {"exp_id": "1"}, FROZEN_CLOCK)
    assert_ok(effects)
    assert effects["store"] is not None
    exp = lookup_experiment(effects["store"], "LM-0001")
    assert exp["status"] == "complete"


def test_complete_does_not_mutate_input():
    """The input store must not be modified."""
    store = copy.deepcopy(MINIMAL_STORE)
    original_status = store["experiments"][0]["status"]
    transform_exp_complete(store, {"exp_id": "LM-0001"}, FROZEN_CLOCK)
    assert store["experiments"][0]["status"] == original_status


# ============================================================
# TESTS: transform_exp_show
# ============================================================

def test_show_minimal():
    """Show experiment with no related records."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_show(store, {"exp_id": "LM-0001"}, FROZEN_CLOCK)
    assert_ok(effects)
    assert effects["store"] is None  # read-only
    text = stdout_text(effects)
    assert "LM-0001" in text
    assert "purification" in text
    assert "active (9d)" in text
    # should NOT have sub-sections
    assert "Strains:" not in text
    assert "Links:" not in text
    assert "Files:" not in text
    assert "Purification details:" not in text
    assert "Assay details:" not in text
    assert "Samples" not in text


def test_show_full():
    """Show experiment with all related records populated."""
    # Need LM-0001 in the store for the link target
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    store["experiments"].append({
        "id": "LM-0001",
        "type": "purification",
        "title": "test-pur",
        "shortname": "test-pur",
        "date_started": "2025-05-15",
        "date_completed": None,
        "status": "active",
        "folder_name": "20250515_LM-0001_pur_test-pur",
        "protocol_id": None,
        "notes": None,
        "created_at": "2025-05-15T00:00:00",
    })
    store["experiments"].append({
        "id": "LM-0003",
        "type": "gelshift",
        "title": "gel-rep",
        "shortname": "gel-rep",
        "date_started": "2025-06-01",
        "date_completed": None,
        "status": "active",
        "folder_name": "20250601_LM-0003_gs_gel-rep",
        "protocol_id": None,
        "notes": None,
        "created_at": "2025-06-01T00:00:00",
    })
    effects = transform_exp_show(store, {"exp_id": "LM-0002"}, FROZEN_CLOCK)
    assert_ok(effects)
    assert effects["store"] is None  # read-only
    text = stdout_text(effects)
    # main fields
    assert "LM-0002" in text
    assert "gelshift" in text
    assert "complete (16d)" in text
    assert "P-001" in text
    assert "[2025-05-20] initial notes" in text
    # sub-sections
    assert "Strains:" in text
    assert "LY456" in text
    assert "(test)" in text
    assert "Links:" in text
    assert "-> LM-0001 (uses_prep_from)" in text
    assert "<- LM-0003 (replicate_of)" in text
    assert "Files:" in text
    assert "gel-01.tiff" in text
    assert "[gel]" in text
    assert "coomassie" in text
    assert "Assay details:" in text
    assert "protein_prep: LM-0001" in text
    assert "dna_template: pUC19" in text
    assert "conditions: 37C 1hr" in text
    assert "Samples (2):" in text
    assert "orc4-mut" in text
    assert "no protein" in text


def test_show_nonexistent():
    """Non-existent experiment -> error effects."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_show(store, {"exp_id": "LM-9999"}, FROZEN_CLOCK)
    assert_fail(effects)
    assert "No experiment found" in effects["stderr"][0]


def test_show_invalid_id():
    """Invalid ID -> error effects."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_show(store, {"exp_id": "garbage"}, FROZEN_CLOCK)
    assert_fail(effects)
    assert "Invalid experiment ID" in effects["stderr"][0]


def test_show_does_not_mutate_store():
    """exp show must not modify the store (it's read-only)."""
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    original = copy.deepcopy(store)
    transform_exp_show(store, {"exp_id": "LM-0002"}, FROZEN_CLOCK)
    assert store == original


# ============================================================
# RUNNER
# ============================================================

def run_all():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    passed = 0
    failed = 0
    for t in tests:
        name = t.__name__
        try:
            t()
            print(f"  PASS  {name}")
            passed += 1
        except Exception as e:
            print(f"  FAIL  {name}: {e}")
            failed += 1
    print()
    print(f"{passed} passed, {failed} failed, {passed + failed} total")
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(run_all())
