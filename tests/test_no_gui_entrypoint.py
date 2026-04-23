import mdpertool.no_gui.no_gui as no_gui_entry


def test_no_gui_main_delegates_to_cli_runner(monkeypatch):
    captured = {}

    def fake_runner(args):
        captured["args"] = args

    monkeypatch.setattr(no_gui_entry, "run_mdpertool_from_cli", fake_runner)
    monkeypatch.setattr(
        "sys.argv",
        [
            "no_gui.py",
            "-p",
            "protein.pdb",
            "-pert_res",
            "ASN17A",
            "-speed_factor",
            "4",
        ],
    )

    exit_code = no_gui_entry.main()

    assert exit_code == 0
    assert captured["args"].topology == "protein.pdb"
    assert captured["args"].perturbed_residues == ["ASN17A"]
    assert captured["args"].velocity_speed_factor == [4]
