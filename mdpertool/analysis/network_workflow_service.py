"""Workflow helpers for analysis network preparation and validation."""

from __future__ import annotations

import os
import csv
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence

from .createRNetwork import MultiTaskEngine
from .pdbsum_conservation_puller import get_conservation_scores


def collect_network_parameters(
    number_of_threads: int,
    pdb_file: str,
    cutoff: float,
    response_time_file: str,
    output_file_name: str,
    output_directory: str,
    source_residue: str,
    node_threshold: int,
    use_node_threshold_condition: bool,
    all_residue_as_target: bool,
) -> Dict[str, Any]:
    """Collect normalized network parameters from UI values."""
    resolved_node_threshold = None if use_node_threshold_condition else node_threshold
    return {
        "number_of_threads": number_of_threads,
        "pdb": str(pdb_file).strip(),
        "cutoff": cutoff,
        "retime_file": str(response_time_file).strip(),
        "outputFileName": output_file_name,
        "output_directory": str(output_directory).strip(),
        "source": source_residue,
        "node_threshold": resolved_node_threshold,
        "node_threshold_use_condition": use_node_threshold_condition,
        "all_residue_as_target": all_residue_as_target,
        "create_output": True,
    }


def _resolve_existing_file_path(path_value: str) -> str:
    """Resolve a file path robustly and ensure it exists.

    Supports absolute and relative paths. Relative paths are resolved against cwd.
    """
    raw = str(path_value or "").strip()
    if not raw:
        raise FileNotFoundError("Required file path is empty.")

    expanded = os.path.expanduser(os.path.expandvars(raw))
    candidate_paths = [expanded]

    if not os.path.isabs(expanded):
        candidate_paths.append(os.path.abspath(expanded))

    for candidate in candidate_paths:
        if os.path.isfile(candidate):
            return os.path.normpath(candidate)

    raise FileNotFoundError(raw)


def validate_network_input_files(network_params: Dict[str, Any]) -> Dict[str, Any]:
    """Validate and normalize topology/response-time file paths before network calculation."""
    pdb_path = _resolve_existing_file_path(network_params.get("pdb", ""))
    response_path = _resolve_existing_file_path(network_params.get("retime_file", ""))

    if os.path.getsize(response_path) == 0:
        raise RuntimeError(
            "Response-time CSV is empty. Please run decomposition/response-time calculation again "
            "or select a valid responseTimes_*.csv file."
        )

    has_data_row = False
    with open(response_path, "r", encoding="utf-8", errors="ignore") as response_file:
        for line in response_file:
            if line.strip():
                has_data_row = True
                break

    if not has_data_row:
        raise RuntimeError(
            "Response-time CSV has no readable data rows. "
            "Please regenerate response times and retry network analysis."
        )

    network_params["pdb"] = pdb_path
    network_params["retime_file"] = response_path
    return network_params


def collect_target_residues(
    all_residue_as_target: bool,
    all_target_items: Sequence[str],
    selected_target_items: Sequence[str],
) -> List[str]:
    """Collect target residues according to target-selection mode."""
    if all_residue_as_target:
        return list(all_target_items)
    return list(selected_target_items)


def read_target_items_from_widgets(
    target_res_combo_box: Any,
    selected_target_residues_list_widget: Any,
) -> Dict[str, List[str]]:
    """Read all/selected target residue labels from UI widgets."""
    all_target_items = [
        target_res_combo_box.itemText(index)[:-1]
        for index in range(target_res_combo_box.count())
    ]
    selected_target_items = [
        selected_target_residues_list_widget.item(index).text()[:-1]
        for index in range(selected_target_residues_list_widget.count())
    ]
    return {
        "all_target_items": all_target_items,
        "selected_target_items": selected_target_items,
    }


def collect_conservation_settings(
    use_conservation: bool,
    pdb_id: str,
    chain_id: str,
    conservation_threshold: float,
    save_conservation_scores: bool = False,
) -> Dict[str, Any]:
    """Collect conservation settings from UI values."""
    return {
        "use_conservation": use_conservation,
        "pdb_id": pdb_id,
        "chain": chain_id,
        "conservation_threshold": conservation_threshold,
        "save_conservation_scores": save_conservation_scores,
    }


def build_network_run_context(
    number_of_threads: int,
    pdb_file: str,
    cutoff: float,
    response_time_file: str,
    output_file_name: str,
    output_directory: str,
    source_residue: str,
    node_threshold: int,
    use_node_threshold_condition: bool,
    all_residue_as_target: bool,
    all_target_items: Sequence[str],
    selected_target_items: Sequence[str],
    use_conservation: bool,
    pdb_id: str,
    chain_id: str,
    conservation_threshold: float,
    atom_pair_checked: bool,
    calpha_checked: bool,
    center_of_mass_checked: bool,
    save_conservation_scores: bool = False,
) -> Dict[str, Any]:
    """Build normalized context for a complete network-calculation run."""
    network_params = collect_network_parameters(
        number_of_threads=number_of_threads,
        pdb_file=pdb_file,
        cutoff=cutoff,
        response_time_file=response_time_file,
        output_file_name=output_file_name,
        output_directory=output_directory,
        source_residue=source_residue,
        node_threshold=node_threshold,
        use_node_threshold_condition=use_node_threshold_condition,
        all_residue_as_target=all_residue_as_target,
    )

    target_residues = collect_target_residues(
        all_residue_as_target=all_residue_as_target,
        all_target_items=all_target_items,
        selected_target_items=selected_target_items,
    )

    conservation_settings = collect_conservation_settings(
        use_conservation=use_conservation,
        pdb_id=pdb_id,
        chain_id=chain_id,
        conservation_threshold=conservation_threshold,
        save_conservation_scores=save_conservation_scores,
    )

    method = determine_network_method(
        atom_pair_checked=atom_pair_checked,
        calpha_checked=calpha_checked,
        center_of_mass_checked=center_of_mass_checked,
    )

    network_output_directory = prepare_network_output_directory(
        output_directory=network_params["output_directory"],
        source_residue=network_params["source"],
    )

    return {
        "network_params": network_params,
        "target_residues": target_residues,
        "conservation_settings": conservation_settings,
        "method": method,
        "network_output_directory": network_output_directory,
    }


def build_network_run_context_from_ui(target: Any) -> Dict[str, Any]:
    """Build network run context directly from UI host object widgets/fields."""
    target_item_payload = read_target_items_from_widgets(
        target_res_combo_box=target.target_res_comboBox,
        selected_target_residues_list_widget=target.selected_target_residues_listWidget,
    )

    return build_network_run_context(
        number_of_threads=target.Number_of_thread_for_network_spinBox.value(),
        pdb_file=target.boundForm_pdb_lineedit.text(),
        cutoff=target.network_cutoff_spinBox.value(),
        response_time_file=target.response_time_lineEdit.text(),
        output_file_name=target.PPI_Network_name_lineedit.text(),
        output_directory=target.net_output_directory_lineedit.text(),
        source_residue=target.source_res_comboBox.currentText()[:-1],
        node_threshold=target.node_threshold_spinBox.value(),
        use_node_threshold_condition=target.node_threshold_checkBox.isChecked(),
        all_residue_as_target=target.all_targets_checkBox.isChecked(),
        all_target_items=target_item_payload["all_target_items"],
        selected_target_items=target_item_payload["selected_target_items"],
        use_conservation=target.use_conservation_checkBox.isChecked(),
        pdb_id=target.conservation_PDB_ID_lineEdit.text(),
        chain_id=target.conservation_pdb_chain_id_lineedit.text(),
        conservation_threshold=target.conserv_score_doubleSpinBox.value(),
        atom_pair_checked=target.atomPair_checkBox.isChecked(),
        calpha_checked=target.Calpha_checkBox.isChecked(),
        center_of_mass_checked=target.center_of_mass_checkBox.isChecked(),
        save_conservation_scores=False,
    )


def apply_network_params(target: Any, network_params: Dict[str, Any]) -> None:
    """Apply normalized network parameters as attributes on target object."""
    for key, value in network_params.items():
        setattr(target, key, value)


def apply_network_run_context(target: Any, run_context: Dict[str, Any]) -> None:
    """Apply full run-context attributes onto host target object."""
    target.network_params = run_context["network_params"]
    target._target_residues_for_calc = run_context["target_residues"]
    target._conservation_settings_for_calc = run_context["conservation_settings"]
    target._network_output_directory = run_context["network_output_directory"]
    apply_network_params(target, target.network_params)


def prepare_general_network_engine_from_ui(target: Any, atom_type: str = "CA") -> MultiTaskEngine:
    """Initialize runtime state, apply UI run context, and build network engine."""
    initialize_network_runtime_state(target)
    run_context = build_network_run_context_from_ui(target)
    apply_network_run_context(target, run_context)

    try:
        target.network_params = validate_network_input_files(target.network_params)
    except FileNotFoundError as file_error:
        missing_path = str(file_error)
        raise RuntimeError(
            "Network preparation failed because an input file could not be found: "
            f"{missing_path}. Please re-select topology and response-time files."
        ) from file_error
    except RuntimeError:
        raise

    return build_network_engine(
        parameters=target.network_params,
        method=run_context["method"],
        output_directory=target._network_output_directory,
        atom_type=atom_type,
    )


def initialize_network_runtime_state(target: Any) -> None:
    """Reset runtime state fields before starting network calculation."""
    target.active_workers = 0
    target._expected_workers = 0
    target._completed_workers = 0
    target._network_finalize_done = False
    target.network_holder = []
    target.log_holder = []
    target.initial_network = None


def prepare_network_ui_for_run(
    target: Any,
    progress_manager: Any,
    label: str = "Calculating network...",
    title: str = "Network Calculation",
) -> Any:
    """Disable UI and open indeterminate progress dialog for network run."""
    target.setEnabled(False)
    return progress_manager.show_indeterminate(
        label=label,
        title=title,
        window_modal=False,
        frameless=True,
        size=(400, 100),
        cancel_button_text=None,
    )


def restore_network_ui_after_run(target: Any, progress_manager: Any) -> None:
    """Close progress UI and re-enable host window after run/error."""
    progress_manager.close(process_events=True)
    target.setEnabled(True)


def start_pair_network_workers(
    engine: Any,
    target_residues: Sequence[str],
    threadpool: Any,
    on_progress: Any,
    on_started: Any,
    on_result: Any,
    on_finished: Any,
) -> int:
    """Start pair-network workers and return expected worker count."""
    engine.run_pair_network_calculation(target_residues)
    expected_workers = len(getattr(engine, "work", []))

    for work in getattr(engine, "work", []):
        work.signals.progress_on_net_calc.connect(on_progress)
        work.signals.work_started.connect(on_started)
        work.signals.result.connect(on_result)
        work.signals.finished.connect(on_finished)
        threadpool.start(work)

    return expected_workers


def wire_general_network_worker(
    thread: Any,
    worker: Any,
    on_ready: Any,
    on_failed: Any,
) -> None:
    """Wire general-network worker signals and lifecycle handlers."""
    thread.started.connect(worker.run)
    worker.finished.connect(on_ready)
    worker.failed.connect(on_failed)
    worker.finished.connect(thread.quit)
    worker.failed.connect(thread.quit)
    worker.finished.connect(worker.deleteLater)
    worker.failed.connect(worker.deleteLater)
    thread.finished.connect(thread.deleteLater)


def start_general_network_background_worker(
    target: Any,
    engine: Any,
    thread_class: Any,
    worker_class: Any,
    on_ready: Any,
    on_failed: Any,
) -> Dict[str, Any]:
    """Create, wire, and start general-network thread/worker; store references on target."""
    thread = thread_class(target)
    worker = worker_class(engine)
    worker.moveToThread(thread)

    wire_general_network_worker(
        thread=thread,
        worker=worker,
        on_ready=on_ready,
        on_failed=on_failed,
    )

    target._general_network_thread = thread
    target._general_network_worker = worker
    thread.start()

    return {
        "thread": thread,
        "worker": worker,
    }


def build_network_engine(
    parameters: Dict[str, Any],
    method: str,
    output_directory: str,
    atom_type: str = "CA",
) -> MultiTaskEngine:
    """Construct and return MultiTaskEngine from normalized parameters."""
    return MultiTaskEngine(
        pdb_file=parameters["pdb"],
        cutoff=parameters["cutoff"],
        re_time_file=parameters["retime_file"],
        source=parameters["source"],
        node_threshold=parameters["node_threshold"],
        output_file_name=parameters["outputFileName"],
        write_outputs=parameters["create_output"],
        output_directory=output_directory,
        method=method,
        atom_type=atom_type,
    )


def determine_network_method(
    atom_pair_checked: bool,
    calpha_checked: bool,
    center_of_mass_checked: bool,
) -> str:
    """Resolve network calculation method from checkbox states."""
    if atom_pair_checked:
        return "any"
    if calpha_checked:
        return "selected_atom"
    if center_of_mass_checked:
        return "center_of_mass"
    return ""


def prepare_network_output_directory(output_directory: str, source_residue: str) -> str:
    """Create and return output directory for current network run."""
    general_output_folder = os.path.join(output_directory, "network_outputs")
    Path(general_output_folder).mkdir(parents=True, exist_ok=True)

    folder_name = f"output_{source_residue}"
    output_folder_directory = os.path.join(general_output_folder, folder_name)
    Path(output_folder_directory).mkdir(parents=True, exist_ok=True)
    return output_folder_directory


def validate_general_network_inputs(
    initial_network: Any,
    residue_ids: Optional[Sequence[str]],
    response_time_count: int,
    source_residue: str,
    target_residues: Sequence[str],
) -> Dict[str, Any]:
    """Validate general network outputs and derive valid/skipped targets."""
    if initial_network is None or residue_ids is None:
        return {
            "status": "network_failed",
            "title": "Network Preparation Failed",
            "message": "General network could not be created. Please verify that the topology and response-time files belong to the same system.",
        }

    if len(residue_ids) != response_time_count:
        return {
            "status": "length_mismatch",
            "title": "Mismatch Error!",
            "message": "The number of residues in the topology file you have provided is not equal to the response time file.",
        }

    available_residues = set(residue_ids)
    if source_residue not in available_residues:
        return {
            "status": "source_missing",
            "title": "Source Residue Not Found",
            "message": f"Selected source residue '{source_residue}' is not present in the loaded topology/response-time dataset.",
        }

    valid_targets = [
        target for target in target_residues if target in available_residues and target != source_residue
    ]
    skipped_targets = [target for target in target_residues if target not in available_residues]

    if not valid_targets:
        return {
            "status": "no_valid_targets",
            "title": "No Valid Target Residues",
            "message": "None of the selected target residues were found in the loaded topology/response-time dataset.",
            "valid_targets": [],
            "skipped_targets": skipped_targets,
        }

    return {
        "status": "ok",
        "valid_targets": valid_targets,
        "skipped_targets": skipped_targets,
    }


def validate_general_network_from_target(
    target: Any,
    initial_network: Any,
    residue_ids: Optional[Sequence[str]],
    response_time_count: int,
) -> Dict[str, Any]:
    """Validate general-network outputs using attributes on UI host object."""
    return validate_general_network_inputs(
        initial_network=initial_network,
        residue_ids=residue_ids,
        response_time_count=response_time_count,
        source_residue=target.source,
        target_residues=target._target_residues_for_calc,
    )


def build_validation_warning_payload(validation_result: Dict[str, Any]) -> Optional[Dict[str, str]]:
    """Build warning title/message payload for non-ok validation results."""
    if validation_result.get("status") == "ok":
        return None

    return {
        "title": str(validation_result.get("title") or "Validation Error"),
        "message": str(validation_result.get("message") or "Invalid network input configuration."),
    }


def build_skipped_targets_warning(skipped_targets: Sequence[str]) -> Optional[str]:
    """Build progress warning text for skipped target residues."""
    if not skipped_targets:
        return None
    return (
        f"Warning: {len(skipped_targets)} target residue(s) were skipped because they are not in the current network."
    )


def evaluate_general_network_ready(validation_result: Dict[str, Any]) -> Dict[str, Any]:
    """Map validation output to UI-ready decision payload for ready phase."""
    warning_payload = build_validation_warning_payload(validation_result)
    if warning_payload is not None:
        return {
            "should_stop": True,
            "warning_payload": warning_payload,
            "valid_targets": [],
            "skipped_targets_warning": None,
        }

    valid_targets = list(validation_result.get("valid_targets", []))
    skipped_targets = list(validation_result.get("skipped_targets", []))
    return {
        "should_stop": False,
        "warning_payload": None,
        "valid_targets": valid_targets,
        "skipped_targets_warning": build_skipped_targets_warning(skipped_targets),
    }


def build_network_error_warning_payload(error_text: str, phase: str = "runtime") -> Dict[str, str]:
    """Build standardized warning payload for network startup/runtime failures."""
    normalized_phase = str(phase).strip().lower()
    if normalized_phase == "startup":
        return {
            "title": "Network Calculation Error",
            "message": f"An error occurred before network calculation started: {error_text}",
        }

    return {
        "title": "Network Calculation Error",
        "message": f"An error occurred during network calculation: {error_text}",
    }


def finalize_network_failure(
    target: Any,
    progress_manager: Any,
    error_text: str,
    phase: str = "runtime",
) -> Dict[str, str]:
    """Restore UI state after failure and return warning payload."""
    restore_network_ui_after_run(target, progress_manager)
    return build_network_error_warning_payload(str(error_text), phase=phase)


def write_conservation_scores(
    output_directory: str,
    residue_ids: Sequence[str],
    conservation_scores: Sequence[float],
    pdb_id: str,
) -> str:
    """Write conservation scores into CSV and return file path."""
    output_path = os.path.join(output_directory, f"conservation_{pdb_id}.csv")
    with open(output_path, "w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        for residue_id, score in zip(residue_ids, conservation_scores):
            writer.writerow((residue_id, score))
    return output_path


def resolve_targets_with_optional_conservation(
    valid_targets: Sequence[str],
    conservation_settings: Dict[str, Any],
    bound_pdb: str,
    output_directory: str,
) -> Dict[str, Any]:
    """Apply optional conservation filter and return target set with metadata."""
    if not conservation_settings.get("use_conservation", False):
        return {
            "targets": list(valid_targets),
            "conservation_used": False,
            "conservation_written": False,
        }

    residue_ids, conservation_scores = get_conservation_scores(
        pdb_id=conservation_settings.get("pdb_id", ""),
        chain_id=conservation_settings.get("chain", ""),
        cutoff=conservation_settings.get("conservation_threshold", 0.0),
        bound_pdb=bound_pdb,
    )

    if conservation_settings.get("save_conservation_scores", False):
        write_conservation_scores(
            output_directory=output_directory,
            residue_ids=residue_ids,
            conservation_scores=conservation_scores,
            pdb_id=conservation_settings.get("pdb_id", "unknown"),
        )

    filtered_targets = list(set(residue_ids).intersection(valid_targets))
    return {
        "targets": filtered_targets,
        "conservation_used": True,
        "conservation_written": bool(conservation_settings.get("save_conservation_scores", False)),
    }


def resolve_targets_and_start_pair_workers(
    engine: Any,
    valid_targets: Sequence[str],
    conservation_settings: Dict[str, Any],
    bound_pdb: str,
    output_directory: str,
    threadpool: Any,
    on_progress: Any,
    on_started: Any,
    on_result: Any,
    on_finished: Any,
) -> Dict[str, Any]:
    """Resolve target residues (optional conservation) and start pair workers."""
    target_resolution = resolve_targets_with_optional_conservation(
        valid_targets=valid_targets,
        conservation_settings=conservation_settings,
        bound_pdb=bound_pdb,
        output_directory=output_directory,
    )

    expected_workers = start_pair_network_workers(
        engine=engine,
        target_residues=target_resolution["targets"],
        threadpool=threadpool,
        on_progress=on_progress,
        on_started=on_started,
        on_result=on_result,
        on_finished=on_finished,
    )

    return {
        "expected_workers": expected_workers,
        "targets": target_resolution["targets"],
        "conservation_used": target_resolution.get("conservation_used", False),
        "conservation_written": target_resolution.get("conservation_written", False),
    }


def apply_pair_worker_start_state(target: Any, engine: Any, worker_start_result: Dict[str, Any]) -> int:
    """Apply pair-worker start state to host object and return expected worker count."""
    target._active_network_engine = engine
    target._completed_workers = 0
    target._network_finalize_done = False

    expected_workers = int(worker_start_result.get("expected_workers", 0))
    target._expected_workers = expected_workers
    return expected_workers


def process_general_network_ready(
    target: Any,
    initial_network: Any,
    residue_ids: Optional[Sequence[str]],
    response_time_count: int,
    engine: Any,
    threadpool: Any,
    on_progress: Any,
    on_started: Any,
    on_result: Any,
    on_finished: Any,
    output_directory: Optional[str] = None,
) -> Dict[str, Any]:
    """Process general-network ready phase and return decision/start payload."""
    validation_result = validate_general_network_from_target(
        target=target,
        initial_network=initial_network,
        residue_ids=residue_ids,
        response_time_count=response_time_count,
    )
    ready_decision = evaluate_general_network_ready(validation_result)
    if ready_decision["should_stop"]:
        return {
            "status": "stop",
            "warning_payload": ready_decision["warning_payload"],
            "skipped_targets_warning": None,
            "worker_start_result": None,
            "expected_workers": 0,
        }

    worker_start_result = resolve_targets_and_start_pair_workers(
        engine=engine,
        valid_targets=ready_decision["valid_targets"],
        conservation_settings=target._conservation_settings_for_calc,
        bound_pdb=target.pdb,
        output_directory=output_directory or getattr(target, "_network_output_directory", target.output_directory),
        threadpool=threadpool,
        on_progress=on_progress,
        on_started=on_started,
        on_result=on_result,
        on_finished=on_finished,
    )

    return {
        "status": "started",
        "warning_payload": None,
        "skipped_targets_warning": ready_decision["skipped_targets_warning"],
        "worker_start_result": worker_start_result,
        "expected_workers": int(worker_start_result.get("expected_workers", 0)),
    }


def finalize_general_network_ready_stop(
    target: Any,
    progress_manager: Any,
    ready_phase_result: Dict[str, Any],
) -> Optional[Dict[str, str]]:
    """Handle ready-phase stop result by restoring UI and returning warning payload."""
    if ready_phase_result.get("status") != "stop":
        return None

    restore_network_ui_after_run(target, progress_manager)
    warning_payload = ready_phase_result.get("warning_payload")
    if isinstance(warning_payload, dict):
        return {
            "title": str(warning_payload.get("title", "Validation Error")),
            "message": str(warning_payload.get("message", "Invalid network input configuration.")),
        }
    return {
        "title": "Validation Error",
        "message": "Invalid network input configuration.",
    }


def apply_general_network_ready_start_state(
    target: Any,
    engine: Any,
    ready_phase_result: Dict[str, Any],
) -> int:
    """Apply worker-start state for a started ready-phase result."""
    if ready_phase_result.get("status") != "started":
        return 0

    worker_start_result = ready_phase_result.get("worker_start_result")
    if not isinstance(worker_start_result, dict):
        return 0

    return apply_pair_worker_start_state(target, engine, worker_start_result)


def finalize_general_network_ready_started(
    target: Any,
    engine: Any,
    progress_manager: Any,
    ready_phase_result: Dict[str, Any],
    on_progress_message: Optional[Any] = None,
) -> int:
    """Finalize started ready-phase by posting warning text and handling no-worker restore."""
    skipped_targets_warning = ready_phase_result.get("skipped_targets_warning")
    if skipped_targets_warning and callable(on_progress_message):
        on_progress_message(skipped_targets_warning)

    expected_workers = apply_general_network_ready_start_state(
        target=target,
        engine=engine,
        ready_phase_result=ready_phase_result,
    )
    if expected_workers == 0:
        restore_network_ui_after_run(target, progress_manager)

    return expected_workers


def handle_general_network_ready_callback(
    target: Any,
    progress_manager: Any,
    initial_network: Any,
    residue_ids: Optional[Sequence[str]],
    response_time_count: int,
    engine: Any,
    threadpool: Any,
    on_progress: Any,
    on_started: Any,
    on_result: Any,
    on_finished: Any,
    on_progress_message: Optional[Any] = None,
    output_directory: Optional[str] = None,
) -> Optional[Dict[str, str]]:
    """Run complete ready-callback flow; return warning payload only when flow should stop."""
    ready_phase_result = process_general_network_ready(
        target=target,
        initial_network=initial_network,
        residue_ids=residue_ids,
        response_time_count=response_time_count,
        engine=engine,
        threadpool=threadpool,
        on_progress=on_progress,
        on_started=on_started,
        on_result=on_result,
        on_finished=on_finished,
        output_directory=output_directory,
    )

    warning_payload = finalize_general_network_ready_stop(
        target=target,
        progress_manager=progress_manager,
        ready_phase_result=ready_phase_result,
    )
    if warning_payload is not None:
        return warning_payload

    finalize_general_network_ready_started(
        target=target,
        engine=engine,
        progress_manager=progress_manager,
        ready_phase_result=ready_phase_result,
        on_progress_message=on_progress_message,
    )
    return None


def present_network_failure_warning(
    target: Any,
    progress_manager: Any,
    error_text: str,
    phase: str,
    show_warning_fn: Any,
    stylesheet: Any,
    warning_presenter: Any,
) -> Dict[str, str]:
    """Finalize network failure and present warning via provided presenter function."""
    warning_payload = finalize_network_failure(
        target=target,
        progress_manager=progress_manager,
        error_text=str(error_text),
        phase=phase,
    )
    warning_presenter(
        show_warning_fn=show_warning_fn,
        owner=target,
        payload=warning_payload,
        stylesheet=stylesheet,
    )
    return warning_payload
