"""Advanced platform option management for OpenMM simulations.

This module provides UI and logic for managing advanced OpenMM platform properties:
- Determinism enforcement (CUDA/OpenCL)
- Advanced options (UseCpuPme, UseBlockingSync)
"""

from PySide2 import QtWidgets, QtCore
from gui.ui_styles import Style


def initialize_advanced_platform_options(main_window):
    """Create and manage advanced platform options in the UI.
    
    Adds programmatically:
    1. Determinism checkbox (for CUDA/OpenCL platforms)
    2. Advanced Options collapsible group
       - UseCpuPme toggle
       - UseBlockingSync toggle
    
    Parameters:
    -----------
    main_window : MainWindow
        Reference to main GUI window
    """
    # Bind widgets that are defined in MAIN_GUI.ui
    main_window.determinism_checkBox = getattr(main_window, "determinism_checkBox", None)
    main_window.use_cpu_pme_checkBox = getattr(main_window, "use_cpu_pme_checkBox", None)
    main_window.use_blocking_sync_checkBox = getattr(main_window, "use_blocking_sync_checkBox", None)
    main_window.optimize_pme_checkBox = getattr(main_window, "optimize_pme_checkBox", None)
    main_window.nonBounded_Method_comboBox = getattr(main_window, "nonBounded_Method_comboBox", None)

    platform_combo = getattr(main_window, "equ_platform_comboBox", None)
    if (not platform_combo or
            main_window.determinism_checkBox is None or
            main_window.use_cpu_pme_checkBox is None or
            main_window.use_blocking_sync_checkBox is None):
        return

    main_window.determinism_checkBox.setToolTip(
        "Enable deterministic force calculations for reproducible simulations.\n"
        "CUDA/OpenCL only. Slightly reduces performance."
    )
    main_window._default_use_cpu_pme_text = main_window.use_cpu_pme_checkBox.text()
    main_window._default_use_cpu_pme_tooltip = (
        "Use CPU for PME (Particle Mesh Ewald) calculations instead of GPU.\n"
        "Useful for CPU-intensive PME or mixed workloads.\n"
        "Available only when Nonbonded Method is PME and platform is CUDA/OpenCL."
    )
    main_window.use_cpu_pme_checkBox.setToolTip(
        main_window._default_use_cpu_pme_tooltip
    )
    if main_window.optimize_pme_checkBox is not None:
        main_window.optimize_pme_checkBox.setToolTip(
            "Automatically benchmark PME cutoff/CPU-PME settings for better performance.\n"
            "Available only when Nonbonded Method is PME and platform is CUDA/OpenCL."
        )
    main_window.use_blocking_sync_checkBox.setToolTip(
        "Allow CPU to sleep while GPU works. Slightly improves perf but blocks other CPU work.\n"
        "CUDA only."
    )
    
    # === CONNECT PLATFORM CHANGES TO UPDATE AVAILABILTY ===
    main_window.equ_platform_comboBox.currentTextChanged.connect(
        lambda _platform_name: _update_advanced_options_availability(main_window)
    )
    if main_window.nonBounded_Method_comboBox is not None:
        main_window.nonBounded_Method_comboBox.currentTextChanged.connect(
            lambda _method_name: _update_advanced_options_availability(main_window)
        )
    if main_window.optimize_pme_checkBox is not None:
        main_window.optimize_pme_checkBox.stateChanged.connect(
            lambda _state: _update_advanced_options_availability(main_window)
        )
    
    # Initial state
    _update_advanced_options_availability(main_window)


def _update_advanced_options_availability(main_window):
    """Update advanced option availability from current platform and nonbonded method.
    
    Parameters:
    -----------
    main_window : MainWindow
        Reference to main GUI window
    Uses:
    - platform combo (CUDA/OpenCL/CPU/Reference)
    - nonbonded method combo (PME/NoCutoff)
    """
    if not hasattr(main_window, 'determinism_checkBox'):
        return

    platform_name = main_window.equ_platform_comboBox.currentText()
    nonbonded_method = None
    if getattr(main_window, 'nonBounded_Method_comboBox', None) is not None:
        nonbonded_method = main_window.nonBounded_Method_comboBox.currentText()
    is_pme_method = nonbonded_method == 'PME'
    
    # Determinism: CUDA only (OpenCL backend rejects this property in our runtime stack).
    is_determinism_available = platform_name == 'CUDA'
    main_window.determinism_checkBox.setEnabled(is_determinism_available)
    if not is_determinism_available:
        main_window.determinism_checkBox.setChecked(False)
    
    # UseCpuPme: CUDA & OpenCL only
    is_cpu_pme_available = platform_name in ['CUDA', 'OpenCL'] and is_pme_method

    optimize_pme_available = (
        main_window.optimize_pme_checkBox is not None and
        platform_name in ['CUDA', 'OpenCL'] and
        is_pme_method
    )
    optimize_pme_active = False
    if main_window.optimize_pme_checkBox is not None:
        main_window.optimize_pme_checkBox.setEnabled(optimize_pme_available)
        if not optimize_pme_available:
            main_window.optimize_pme_checkBox.setChecked(False)
        optimize_pme_active = main_window.optimize_pme_checkBox.isChecked()

    main_window.use_cpu_pme_checkBox.setEnabled(is_cpu_pme_available and not optimize_pme_active)
    if not is_cpu_pme_available or optimize_pme_active:
        main_window.use_cpu_pme_checkBox.setChecked(False)

    if optimize_pme_active:
        main_window.use_cpu_pme_checkBox.setText("Use CPU PME (managed by Optimize PME)")
        main_window.use_cpu_pme_checkBox.setToolTip(
            "Disabled while Optimize PME is enabled.\n"
            "Optimize PME benchmarks CPU/GPU PME paths and picks the faster one automatically.\n"
            "On OpenCL, DeterministicForces is intentionally not set."
        )
    else:
        default_text = getattr(main_window, '_default_use_cpu_pme_text', 'Use CPU PME')
        default_tooltip = getattr(
            main_window,
            '_default_use_cpu_pme_tooltip',
            "Use CPU for PME (Particle Mesh Ewald) calculations instead of GPU.",
        )
        main_window.use_cpu_pme_checkBox.setText(default_text)
        main_window.use_cpu_pme_checkBox.setToolTip(default_tooltip)
    
    # UseBlockingSync: CUDA only
    is_blocking_sync_available = platform_name == 'CUDA'
    main_window.use_blocking_sync_checkBox.setEnabled(is_blocking_sync_available)
    if not is_blocking_sync_available:
        main_window.use_blocking_sync_checkBox.setChecked(False)


def get_advanced_platform_properties(main_window):
    """Return a dict of advanced platform properties based on UI selections.
    
    Returns:
    --------
    dict : Platform-specific advanced properties for OpenMM
        Keys: 'DeterministicForces', 'UseCpuPme', 'UseBlockingSync'
        Values: 'true' or 'false'
    
    Example:
    --------
    >>> props = get_advanced_platform_properties(main_window)
    >>> # {'DeterministicForces': 'true', 'UseCpuPme': 'false', 'UseBlockingSync': 'false'}
    """
    if not hasattr(main_window, 'determinism_checkBox'):
        return {}
    
    properties = {}
    
    # DeterministicForces (CUDA only)
    if main_window.determinism_checkBox.isEnabled() and main_window.determinism_checkBox.isChecked():
        properties['DeterministicForces'] = 'true'
    
    # UseCpuPme (CUDA/OpenCL)
    if main_window.use_cpu_pme_checkBox.isEnabled() and main_window.use_cpu_pme_checkBox.isChecked():
        properties['UseCpuPme'] = 'true'
    
    # UseBlockingSync (CUDA only)
    if main_window.use_blocking_sync_checkBox.isEnabled() and main_window.use_blocking_sync_checkBox.isChecked():
        properties['UseBlockingSync'] = 'true'
    
    return properties
