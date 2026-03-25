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
    main_window.use_cpu_pme_checkBox.setToolTip(
        "Use CPU for PME (Particle Mesh Ewald) calculations instead of GPU.\n"
        "Useful for CPU-intensive PME or mixed workloads."
    )
    main_window.use_blocking_sync_checkBox.setToolTip(
        "Allow CPU to sleep while GPU works. Slightly improves perf but blocks other CPU work.\n"
        "CUDA only."
    )
    
    # === CONNECT PLATFORM CHANGES TO UPDATE AVAILABILTY ===
    main_window.equ_platform_comboBox.currentTextChanged.connect(
        lambda platform_name: _update_advanced_options_availability(main_window, platform_name)
    )
    
    # Initial state
    _update_advanced_options_availability(main_window, main_window.equ_platform_comboBox.currentText())


def _update_advanced_options_availability(main_window, platform_name):
    """Update which advanced options are enabled based on selected platform.
    
    Parameters:
    -----------
    main_window : MainWindow
        Reference to main GUI window
    platform_name : str
        Name of selected OpenMM platform (CUDA, OpenCL, CPU, Reference)
    """
    if not hasattr(main_window, 'determinism_checkBox'):
        return
    
    # Determinism: CUDA & OpenCL only
    is_determinism_available = platform_name in ['CUDA', 'OpenCL']
    main_window.determinism_checkBox.setEnabled(is_determinism_available)
    if not is_determinism_available:
        main_window.determinism_checkBox.setChecked(False)
    
    # UseCpuPme: CUDA & OpenCL only
    is_cpu_pme_available = platform_name in ['CUDA', 'OpenCL']
    main_window.use_cpu_pme_checkBox.setEnabled(is_cpu_pme_available)
    if not is_cpu_pme_available:
        main_window.use_cpu_pme_checkBox.setChecked(False)
    
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
    
    # DeterministicForces (CUDA/OpenCL)
    if main_window.determinism_checkBox.isEnabled() and main_window.determinism_checkBox.isChecked():
        properties['DeterministicForces'] = 'true'
    
    # UseCpuPme (CUDA/OpenCL)
    if main_window.use_cpu_pme_checkBox.isEnabled() and main_window.use_cpu_pme_checkBox.isChecked():
        properties['UseCpuPme'] = 'true'
    
    # UseBlockingSync (CUDA only)
    if main_window.use_blocking_sync_checkBox.isEnabled() and main_window.use_blocking_sync_checkBox.isChecked():
        properties['UseBlockingSync'] = 'true'
    
    return properties
