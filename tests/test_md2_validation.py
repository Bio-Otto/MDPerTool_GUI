import pytest

from mdpertool.no_gui.MD_2 import _validate_reference_inputs


BASE_ARGS = {
    "pdb_path": "protein.pdb",
    "state_file": "state.xml",
    "protein_ff": "amber96",
    "water_ff": "tip3p",
    "time_step": 1,
    "nonbondedCutoff": 12.0,
    "use_switching_distance": False,
    "switching_distance": None,
    "reference_total_Steps": 1000,
    "platform_name": "CPU",
    "output_directory": "output",
}


def test_validate_reference_inputs_accepts_minimum_valid_args():
    _validate_reference_inputs(**BASE_ARGS)


def test_validate_reference_inputs_rejects_missing_required_fields():
    args = dict(BASE_ARGS)
    args["protein_ff"] = None
    args["output_directory"] = None

    with pytest.raises(ValueError) as exc:
        _validate_reference_inputs(**args)

    assert "protein_ff" in str(exc.value)
    assert "output_directory" in str(exc.value)


def test_validate_reference_inputs_requires_switching_distance_when_enabled():
    args = dict(BASE_ARGS)
    args["use_switching_distance"] = True
    args["switching_distance"] = None

    with pytest.raises(ValueError) as exc:
        _validate_reference_inputs(**args)

    assert "switching_distance" in str(exc.value)
