"""UI Layout Manager - Handle panel visibility and layout controls."""


class UILayoutManager:
    """Manage visibility and toggling of UI analysis panel layouts."""

    def show_visualization_settings_on_analysis(self, analysis_settings_groupBox, show_navigation_button,
                                                hide_navigation_button):
        """Show the visualization settings panel on analysis tab."""
        analysis_settings_groupBox.show()
        show_navigation_button.hide()
        hide_navigation_button.show()

    def hide_visualization_settings_on_analysis(self, show_navigation_button, hide_navigation_button,
                                                analysis_settings_groupBox):
        """Hide the visualization settings panel on analysis tab."""
        analysis_settings_groupBox.hide()
        show_navigation_button.show()
        hide_navigation_button.hide()

    def visualization_Handle_buttons_changing_on_analysis(self, analysis_settings_groupBox, show_navigation_button,
                                                          hide_navigation_button):
        """Trigger visualization settings panel toggle."""
        self.hide_visualization_settings_on_analysis(analysis_settings_groupBox,
                                                     show_navigation_button,
                                                     hide_navigation_button)

    def Handle_Buttons_on_analysis(self, analysis_settings_groupBox, show_navigation_button, hide_navigation_button):
        """Connect show/hide navigation buttons to visualization panel toggle."""
        show_navigation_button.clicked.connect(lambda:
                                               self.show_visualization_settings_on_analysis(
                                                   analysis_settings_groupBox,
                                                   show_navigation_button,
                                                   hide_navigation_button))
        hide_navigation_button.clicked.connect(lambda:
                                               self.hide_visualization_settings_on_analysis(
                                                   show_navigation_button,
                                                   hide_navigation_button,
                                                   analysis_settings_groupBox))

    def show_figure_options_on_analysis(self, figure_settings_groupBox_on_analysis):
        """Show figure export options panel."""
        figure_settings_groupBox_on_analysis.show()

    def hide_figure_options_on_analysis(self, figure_settings_groupBox_on_analysis):
        """Hide figure export options panel."""
        figure_settings_groupBox_on_analysis.hide()

    def Handle_Save_Figure_Options_on_analysis_Changed(self, figure_settings_groupBox_on_analysis):
        """Handle figure options panel visibility changes."""
        self.hide_figure_options_on_analysis(figure_settings_groupBox_on_analysis)

    def Handle_Save_Figure_Options_on_analysis(self, save_as_png_pushButton, hide_figure_settings_pushButton,
                                               width_horizontalSlider, height_horizontalSlider, dpi_horizontalSlider,
                                               ray_horizontalSlider, figure_settings_groupBox_on_analysis,
                                               pymol_width_label, pymol_height_label, pymol_dpi_label, pymol_ray_label):
        """Connect figure export option controls to their respective handlers."""
        save_as_png_pushButton.clicked.connect(
            lambda: self.show_figure_options_on_analysis(figure_settings_groupBox_on_analysis))
        hide_figure_settings_pushButton.clicked.connect(
            lambda: self.hide_figure_options_on_analysis(figure_settings_groupBox_on_analysis))
        width_horizontalSlider.valueChanged.connect(
            lambda: self.figure_width_label_on_analysis(width_horizontalSlider, pymol_width_label))
        height_horizontalSlider.valueChanged.connect(
            lambda: self.figure_height_label_on_analysis(height_horizontalSlider, pymol_height_label))
        dpi_horizontalSlider.valueChanged.connect(
            lambda: self.figure_dpi_label_on_analysis(dpi_horizontalSlider, pymol_dpi_label))
        ray_horizontalSlider.valueChanged.connect(
            lambda: self.figure_ray_label_on_analysis(ray_horizontalSlider, pymol_ray_label))

    def figure_width_label_on_analysis(self, width_horizontalSlider, with_label):
        """Update figure width label when slider changes."""
        with_label.setText("Width: " + str(width_horizontalSlider.value()))

    def figure_height_label_on_analysis(self, height_horizontalSlider, height_label):
        """Update figure height label when slider changes."""
        height_label.setText("Height: " + str(height_horizontalSlider.value()))

    def figure_dpi_label_on_analysis(self, dpi_horizontalSlider, dpi_label):
        """Update figure DPI label when slider changes."""
        dpi_label.setText("Dpi: " + str(dpi_horizontalSlider.value()))

    def figure_ray_label_on_analysis(self, ray_horizontalSlider, ray_label):
        """Update figure ray label when slider changes."""
        ray_label.setText("Ray: " + str(ray_horizontalSlider.value()))
