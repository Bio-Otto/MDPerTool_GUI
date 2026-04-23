from pathlib import Path

from mdpertool.no_gui.write_outputs import write_folder, lof_file_settings


def test_write_folder_creates_directory(tmp_path):
    output_root = tmp_path / "output"
    output_root.mkdir()

    created_path, folder_name = write_folder(
        directory=str(output_root),
        top="protein.pdb",
        res="SER247",
        ff="amber96_tip3p",
        just_min_or_md=False,
    )

    assert folder_name.startswith("protein_SER247_amber96_tip3p_MD_")
    assert Path(created_path).exists()
    assert Path(created_path).is_dir()


def test_log_file_settings_disables_propagation(tmp_path):
    config = lof_file_settings(str(tmp_path))

    assert config["disable_existing_loggers"] is False
    assert config["loggers"]["api_logger"]["propagate"] is False
    assert config["loggers"]["batch_process_logger"]["propagate"] is False
