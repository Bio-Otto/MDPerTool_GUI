"""PyMOL Visualization Manager - Handle 3D structure visualization."""

from PySide2.QtWidgets import QFileDialog


class PyMOLVisualizer:
    """Manage PyMOL 3D visualization and figure export operations."""

    def activate_navigation_on_Pymol(self, created_pymol_widget):
        """Activate navigation tool in PyMOL widget."""
        created_pymol_widget.activate_navigation_tool()
        created_pymol_widget.paintGL()
        created_pymol_widget.update()
        created_pymol_widget.show()

    def deactivate_navigation_on_Pymol(self, created_pymol_widget):
        """Deactivate navigation tool in PyMOL widget."""
        created_pymol_widget.deactivate_navigation_tool()
        created_pymol_widget.paintGL()
        created_pymol_widget.update()
        created_pymol_widget.show()

    def clear_residue_labels(self, created_pymol_widget):
        """Clear all residue labels from PyMOL display."""
        created_pymol_widget.clear_all_labels()
        created_pymol_widget.update()

    def show_beautiful_in_Pymol(self, created_pymol_widget):
        """Display secondary structure in PyMOL."""
        created_pymol_widget.set_ss_figure()
        created_pymol_widget.update()
        created_pymol_widget.show()

    def save_as_png_Pymol(self, main_window, created_pymol_widget, width_horizontalSlider, height_horizontalSlider,
                          dpi_horizontalSlider, ray_horizontalSlider):
        """Export PyMOL visualization as PNG file."""
        # Lazy import to avoid circular dependencies
        from ..message import Message_Boxes
        from ..gui.ui_styles import Style
        
        filedialog = QFileDialog(main_window)
        filedialog.setDefaultSuffix("png")
        filedialog.setNameFilter("PNG Files (*.png);;All files (*.*)")
        filedialog.setAcceptMode(QFileDialog.AcceptSave)
        selected = filedialog.exec()

        if selected:
            filename = filedialog.selectedFiles()[0]
        else:
            return
        
        if filename == "":
            Message_Boxes.Warning_message(main_window, 'png save failed!', "No file name selected.",
                                          Style.MessageBox_stylesheet)
            return

        try:
            created_pymol_widget.get_png_figure(filename, width=width_horizontalSlider.value(),
                                                height=height_horizontalSlider.value(),
                                                dpi=dpi_horizontalSlider.value(),
                                                ray=ray_horizontalSlider.value())
            created_pymol_widget.update()

        except Exception as save_err:
            Message_Boxes.Critical_message(main_window, 'png save failed!', str(save_err), Style.MessageBox_stylesheet)
