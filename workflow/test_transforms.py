"""Golden tests for Phase 2 and Phase 3A transforms.

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
    transform_exp_update,
    transform_exp_link,
    transform_exp_addstrain,
    transform_exp_find,
    transform_exp_list,
    transform_strain_add,
    transform_strain_update,
    transform_strain_show,
    transform_strain_list,
    transform_stage_list,
    transform_status,
    effects_ok,
    effects_fail,
    lookup_experiment,
    lookup_strain,
    append_note,
    format_size,
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
# TESTS: append_note (shared helper)
# ============================================================

def test_append_note_first():
    result = append_note(None, "first entry", FROZEN_CLOCK)
    assert result == "[2025-06-10] first entry"


def test_append_note_append():
    result = append_note("[2025-06-01] old", "new entry", FROZEN_CLOCK)
    assert result == "[2025-06-01] old\n[2025-06-10] new entry"


# ============================================================
# TESTS: transform_exp_update
# ============================================================

def test_update_simple_field():
    """Update a simple experiments-table field (title)."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_update(
        store, {"exp_id": "LM-0001", "field": "title", "value": "new-title"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    exp = lookup_experiment(effects["store"], "LM-0001")
    assert exp["title"] == "new-title"
    text = stdout_text(effects)
    assert "Updated LM-0001" in text
    assert "title:" in text


def test_update_notes_append():
    """Notes field appends with timestamp."""
    store = copy.deepcopy(MINIMAL_STORE)
    store["experiments"][0]["notes"] = "[2025-06-01] old note"
    effects = transform_exp_update(
        store, {"exp_id": "LM-0001", "field": "notes", "value": "new note"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    exp = lookup_experiment(effects["store"], "LM-0001")
    assert "[2025-06-10] new note" in exp["notes"]
    assert "[2025-06-01] old note" in exp["notes"]


def test_update_protected_field():
    """Protected fields are rejected."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_update(
        store, {"exp_id": "LM-0001", "field": "status", "value": "complete"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "cannot be modified" in effects["stderr"][0]


def test_update_unknown_field():
    """Unknown fields are rejected with a list of valid fields."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_update(
        store, {"exp_id": "LM-0001", "field": "bogus", "value": "x"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "Unknown field" in effects["stderr"][0]


def test_update_protein_creates_purification():
    """Setting protein on exp with no purification record creates one."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_update(
        store, {"exp_id": "LM-0001", "field": "protein", "value": "Orc4"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    pur = [r for r in effects["store"]["purification_details"] if r["experiment_id"] == "LM-0001"]
    assert len(pur) == 1
    assert pur[0]["protein"] == "Orc4"
    text = stdout_text(effects)
    assert "Created purification record" in text


def test_update_non_protein_pur_field_fails_without_record():
    """Setting a non-protein purification field without existing record fails."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_update(
        store, {"exp_id": "LM-0001", "field": "yield_mg", "value": "1.5"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "Set 'protein' first" in effects["stderr"][1]


def test_update_assay_field_creates_skeleton():
    """Setting an assay field on exp with no assay record creates one."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_update(
        store, {"exp_id": "LM-0001", "field": "dna_template", "value": "pUC19"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    assay = [r for r in effects["store"]["assay_details"] if r["experiment_id"] == "LM-0001"]
    assert len(assay) == 1
    assert assay[0]["dna_template"] == "pUC19"


def test_update_nonexistent_exp():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_update(
        store, {"exp_id": "LM-9999", "field": "title", "value": "x"}, FROZEN_CLOCK
    )
    assert_fail(effects)


# ============================================================
# TESTS: transform_exp_link
# ============================================================

def _two_exp_store():
    """Store with two experiments for link testing."""
    store = copy.deepcopy(MINIMAL_STORE)
    store["experiments"].append({
        "id": "LM-0002",
        "type": "gelshift",
        "title": "gel-test",
        "shortname": "gel-test",
        "date_started": "2025-06-05",
        "date_completed": None,
        "status": "active",
        "folder_name": "20250605_LM-0002_gs_gel-test",
        "protocol_id": None,
        "notes": None,
        "created_at": "2025-06-05T00:00:00",
    })
    return store


def test_link_simple():
    store = _two_exp_store()
    effects = transform_exp_link(
        store, {"from_id": "LM-0002", "to_id": "LM-0001", "relationship": "follow_up_to"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    links = effects["store"]["experiment_links"]
    assert len(links) == 1
    assert links[0]["from_experiment"] == "LM-0002"
    assert links[0]["to_experiment"] == "LM-0001"


def test_link_uses_prep_from_auto_populates_assay():
    """uses_prep_from creates assay_details with protein_prep."""
    store = _two_exp_store()
    effects = transform_exp_link(
        store, {"from_id": "LM-0002", "to_id": "LM-0001", "relationship": "uses_prep_from"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    assay = [r for r in effects["store"]["assay_details"] if r["experiment_id"] == "LM-0002"]
    assert len(assay) == 1
    assert assay[0]["protein_prep"] == "LM-0001"


def test_link_uses_prep_from_existing_assay_null_prep():
    """uses_prep_from updates existing assay if protein_prep is NULL."""
    store = _two_exp_store()
    store["assay_details"].append({
        "experiment_id": "LM-0002",
        "protein_prep": None,
        "dna_template": "pUC19",
        "dna_prep_method": None,
        "conditions": None,
    })
    effects = transform_exp_link(
        store, {"from_id": "LM-0002", "to_id": "LM-0001", "relationship": "uses_prep_from"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    assay = [r for r in effects["store"]["assay_details"] if r["experiment_id"] == "LM-0002"]
    assert len(assay) == 1
    assert assay[0]["protein_prep"] == "LM-0001"
    assert assay[0]["dna_template"] == "pUC19"  # preserved


def test_link_uses_prep_from_existing_assay_with_prep():
    """uses_prep_from does NOT overwrite existing non-NULL protein_prep."""
    store = _two_exp_store()
    store["assay_details"].append({
        "experiment_id": "LM-0002",
        "protein_prep": "LM-0099",
        "dna_template": None,
        "dna_prep_method": None,
        "conditions": None,
    })
    effects = transform_exp_link(
        store, {"from_id": "LM-0002", "to_id": "LM-0001", "relationship": "uses_prep_from"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    assay = [r for r in effects["store"]["assay_details"] if r["experiment_id"] == "LM-0002"]
    assert assay[0]["protein_prep"] == "LM-0099"  # NOT overwritten


def test_link_duplicate():
    store = _two_exp_store()
    store["experiment_links"].append({
        "from_experiment": "LM-0002",
        "to_experiment": "LM-0001",
        "relationship": "follow_up_to",
    })
    effects = transform_exp_link(
        store, {"from_id": "LM-0002", "to_id": "LM-0001", "relationship": "follow_up_to"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "already exists" in effects["stderr"][0]


def test_link_unknown_relationship_warns():
    """Unknown relationship produces a warning but succeeds."""
    store = _two_exp_store()
    effects = transform_exp_link(
        store, {"from_id": "LM-0002", "to_id": "LM-0001", "relationship": "custom_rel"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    text = stdout_text(effects)
    assert "not a standard relationship" in text
    assert "Linked:" in text


# ============================================================
# TESTS: transform_exp_addstrain
# ============================================================

def test_addstrain_happy():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    # LM-0002 and LY456 exist; add the link
    # First remove existing link
    store["experiment_strains"] = []
    effects = transform_exp_addstrain(
        store, {"exp_id": "LM-0002", "strain_id": "LY456", "role": "test"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    es = effects["store"]["experiment_strains"]
    assert len(es) == 1
    assert es[0]["strain_id"] == "LY456"
    assert es[0]["role"] == "test"


def test_addstrain_duplicate():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_exp_addstrain(
        store, {"exp_id": "LM-0002", "strain_id": "LY456"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "already linked" in effects["stderr"][0]


def test_addstrain_missing_strain():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_addstrain(
        store, {"exp_id": "LM-0001", "strain_id": "NOSUCH"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "No strain found" in effects["stderr"][0]


# ============================================================
# TESTS: transform_strain_add
# ============================================================

def test_strain_add_happy():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_strain_add(
        store, {"strain_id": "LY789", "genotype": "MATa wild-type"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    strain = lookup_strain(effects["store"], "LY789")
    assert strain is not None
    assert strain["genotype"] == "MATa wild-type"
    text = stdout_text(effects)
    assert "Registered: LY789" in text


def test_strain_add_duplicate():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_strain_add(
        store, {"strain_id": "LY456", "genotype": "dup"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "already exists" in effects["stderr"][0]


def test_strain_add_invalid_name():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_strain_add(
        store, {"strain_id": "-bad", "genotype": "x"}, FROZEN_CLOCK
    )
    assert_fail(effects)


def test_strain_add_atypical_format_warns():
    """Non-standard ID format produces a note but succeeds."""
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_strain_add(
        store, {"strain_id": "weird-name", "genotype": "x"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    text = stdout_text(effects)
    assert "doesn't follow typical" in text


# ============================================================
# TESTS: transform_strain_update
# ============================================================

def test_strain_update_label():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_strain_update(
        store, {"strain_id": "LY456", "field": "label", "value": "new-label"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    strain = lookup_strain(effects["store"], "LY456")
    assert strain["label"] == "new-label"
    text = stdout_text(effects)
    assert "Updated LY456" in text


def test_strain_update_notes_append():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    store["strains"][0]["notes"] = "[2025-06-01] old"
    effects = transform_strain_update(
        store, {"strain_id": "LY456", "field": "notes", "value": "new"}, FROZEN_CLOCK
    )
    assert_ok(effects)
    strain = lookup_strain(effects["store"], "LY456")
    assert "[2025-06-10] new" in strain["notes"]
    assert "[2025-06-01] old" in strain["notes"]
    text = stdout_text(effects)
    assert "appended entry" in text


def test_strain_update_unknown_field():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_strain_update(
        store, {"strain_id": "LY456", "field": "bogus", "value": "x"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "Unknown field" in effects["stderr"][0]


def test_strain_update_nonexistent():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_strain_update(
        store, {"strain_id": "NOSUCH", "field": "label", "value": "x"}, FROZEN_CLOCK
    )
    assert_fail(effects)
    assert "No strain found" in effects["stderr"][0]


# ============================================================
# TESTS: transform_strain_show
# ============================================================

def test_strain_show_happy():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_strain_show(store, {"strain_id": "LY456"}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "LY456" in text
    assert "orc4-mut" in text
    assert "MATa orc4-R267A" in text
    assert "KanMX" in text


def test_strain_show_with_experiments():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_strain_show(store, {"strain_id": "LY456"}, FROZEN_CLOCK)
    text = stdout_text(effects)
    assert "Experiments:" in text
    assert "LM-0002" in text
    assert "(test)" in text


def test_strain_show_nonexistent():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_strain_show(store, {"strain_id": "NOSUCH"}, FROZEN_CLOCK)
    assert_fail(effects)


# ============================================================
# TESTS: transform_strain_list
# ============================================================

def test_strain_list_happy():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_strain_list(store, {"no_label_only": False}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "LY456" in text
    assert "orc4-mut" in text
    assert "1 strain" in text


def test_strain_list_empty():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_strain_list(store, {"no_label_only": False}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "No strains registered" in text


def test_strain_list_no_label_filter():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_strain_list(store, {"no_label_only": True}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "All strains have labels" in text


def test_strain_list_no_label_found():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    store["strains"][0]["label"] = None
    effects = transform_strain_list(store, {"no_label_only": True}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "*NO LABEL*" in text
    assert "without labels" in text


# ============================================================
# TESTS: transform_exp_find
# ============================================================

def test_find_all():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_find(store, {"pos_args": []}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "LM-0001" in text
    assert "1 experiment" in text


def test_find_by_query():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_find(store, {"pos_args": ["test-pur"]}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "LM-0001" in text


def test_find_case_insensitive():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_find(store, {"pos_args": ["TEST-PUR"]}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "LM-0001" in text


def test_find_no_results():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_find(store, {"pos_args": ["nonexistent"]}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "No experiments found" in text


def test_find_by_type():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_find(store, {"pos_args": ["--type", "purification"]}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "LM-0001" in text


def test_find_by_type_no_match():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_find(store, {"pos_args": ["--type", "gelshift"]}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "No experiments found" in text


def test_find_invalid_type():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_find(store, {"pos_args": ["--type", "bogus"]}, FROZEN_CLOCK)
    assert_fail(effects)
    assert "Unknown type" in effects["stderr"][0]


def test_find_by_strain():
    store = copy.deepcopy(COMPLETE_EXPERIMENT_STORE)
    effects = transform_exp_find(store, {"pos_args": ["--strain", "LY456"]}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "LM-0002" in text


def test_find_by_nonexistent_strain():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_find(store, {"pos_args": ["--strain", "NOSUCH"]}, FROZEN_CLOCK)
    assert_fail(effects)
    assert "not found in database" in effects["stderr"][0]


# ============================================================
# TESTS: transform_exp_list
# ============================================================

def test_list_all():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_list(store, {"pos_args": []}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "purification (1)" in text
    assert "LM-0001" in text
    assert "1 experiment" in text


def test_list_empty():
    store = copy.deepcopy(MINIMAL_STORE)
    store["experiments"] = []
    effects = transform_exp_list(store, {"pos_args": []}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "No experiments registered" in text


def test_list_with_purification_detail():
    store = copy.deepcopy(MINIMAL_STORE)
    store["purification_details"].append({
        "experiment_id": "LM-0001",
        "protein": "Orc4",
        "induction_od": None, "galactose_hours": None, "columns": None,
        "v5_depletion": 0, "fractions_pooled": None, "yield_mg": 2.5, "storage": None,
    })
    effects = transform_exp_list(store, {"pos_args": []}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "Orc4" in text
    assert "2.5 mg" in text


def test_list_by_type():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_list(store, {"pos_args": ["--type", "purification"]}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "LM-0001" in text


def test_list_invalid_type():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_exp_list(store, {"pos_args": ["--type", "bogus"]}, FROZEN_CLOCK)
    assert_fail(effects)


# ============================================================
# TESTS: transform_stage_list
# ============================================================

def test_stage_list_dir_missing():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_stage_list(store, {"file_list": None}, FROZEN_CLOCK)
    assert_fail(effects)
    assert "does not exist" in effects["stderr"][0]


def test_stage_list_empty():
    store = copy.deepcopy(MINIMAL_STORE)
    effects = transform_stage_list(store, {"file_list": []}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "Staging is empty" in text


def test_stage_list_with_files():
    store = copy.deepcopy(MINIMAL_STORE)
    file_list = [
        {"name": "lemr LM-0001 s1 nanodrop.csv", "size": 2048},
        {"name": "random.txt", "size": 512},
    ]
    effects = transform_stage_list(store, {"file_list": file_list}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "2 files" in text
    assert "LM-0001" in text
    assert "random.txt" in text


def test_stage_list_with_pending():
    store = copy.deepcopy(MINIMAL_STORE)
    store["stage_expectations"].append({
        "id": 1, "experiment_id": "LM-0001", "slot": 1, "descriptor": "gel-01",
        "status": "pending", "created_at": "2025-06-10", "filed_at": None,
    })
    file_list = [{"name": "lemr LM-0001 s1 nanodrop.csv", "size": 1024}]
    effects = transform_stage_list(store, {"file_list": file_list}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "Pending expectations:" in text
    assert "(matched)" in text


# ============================================================
# TESTS: transform_status
# ============================================================

def test_status_minimal():
    store = copy.deepcopy(MINIMAL_STORE)
    fs = {"staging_files": [], "staging_total_size": 0, "lab_notes_date": None}
    effects = transform_status(store, {"fs": fs}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "Experiments:" in text
    assert "1 total" in text
    assert "1 active" in text
    assert "Staging:" in text
    assert "clean" in text
    assert "Strains:" in text
    assert "Samples:" in text


def test_status_with_warnings():
    store = copy.deepcopy(MINIMAL_STORE)
    fs = {"staging_files": ["a.txt"], "staging_total_size": 100, "lab_notes_date": "2025-06-09"}
    effects = transform_status(store, {"fs": fs}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "Attention needed:" in text
    assert "no strains linked" in text
    assert "no notes" in text
    assert "1 file" in text
    assert "Lab notes last modified: 2025-06-09" in text


def test_status_empty():
    store = copy.deepcopy(MINIMAL_STORE)
    store["experiments"] = []
    fs = {"staging_files": [], "staging_total_size": 0, "lab_notes_date": None}
    effects = transform_status(store, {"fs": fs}, FROZEN_CLOCK)
    assert_ok(effects)
    text = stdout_text(effects)
    assert "0 total" in text
    assert "No experiments yet" in text


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
