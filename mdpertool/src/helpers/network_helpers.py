"""Network Analysis Parameters Manager - Handle OpenMM platform and parameter configuration."""

import openmm


class NetworkParametersManager:
    """Manage network analysis parameters and platform detection."""

    def available_platforms(self):
        """Detect and return available OpenMM platforms on system.
        
        Returns:
            Tuple of (platform_names_keys, platform_speed_dict)
        """
        avail_plt_and_speed = dict()

        for index in range(openmm.Platform.getNumPlatforms()):
            platform_name = openmm.Platform.getPlatform(index).getName()
            platform_speed = openmm.Platform.getPlatform(index).getSpeed()
            avail_plt_and_speed[platform_name] = platform_speed

        return avail_plt_and_speed.keys(), avail_plt_and_speed

    def initialize_parameters(self, Number_of_thread_for_network_spinBox, boundForm_pdb_lineedit,
                             network_cutoff_spinBox, response_time_lineEdit, PPI_Network_name_lineedit,
                             net_output_directory_lineedit, source_res_comboBox, node_threshold_spinBox,
                             node_threshold_checkBox, all_targets_checkBox):
        """Initialize network analysis parameters from UI controls.
        
        Args:
            Various Qt widgets containing user-configured parameters
        
        Returns:
            List of (parameter_name, parameter_value) tuples
        """
        number_of_threads = Number_of_thread_for_network_spinBox.value()
        pdb = boundForm_pdb_lineedit.text()
        cutoff = network_cutoff_spinBox.value()
        retime_file = response_time_lineEdit.text()
        outputFileName = PPI_Network_name_lineedit.text()
        output_directory = net_output_directory_lineedit.text()
        source = source_res_comboBox.currentText()[:-1]
        node_threshold = node_threshold_spinBox.value()
        node_threshold_use_condition = node_threshold_checkBox.isChecked()
        all_residue_as_target = all_targets_checkBox.isChecked()
        create_output = True

        if node_threshold_use_condition:
            node_threshold = None

        return [
            ('number_of_threads', number_of_threads),
            ('pdb', pdb),
            ('cutoff', cutoff),
            ('retime_file', retime_file),
            ('outputFileName', outputFileName),
            ('output_directory', output_directory),
            ('source', source),
            ('node_threshold', node_threshold),
            ('node_threshold_use_condition', node_threshold_use_condition),
            ('all_residue_as_target', all_residue_as_target),
            ('create_output', create_output)
        ]
