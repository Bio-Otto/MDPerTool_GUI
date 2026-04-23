"""Progress dialog management helpers for reusable UI feedback flows."""

from __future__ import annotations

from typing import Optional, Tuple

from PySide2 import QtWidgets
from PySide2.QtCore import Qt


class ProgressDialogManager:
    """Create and manage a single reusable QProgressDialog instance."""

    def __init__(self, parent=None, stylesheet: Optional[str] = None):
        self.parent = parent
        self.stylesheet = stylesheet
        self.dialog: Optional[QtWidgets.QProgressDialog] = None

    def show(
        self,
        label: str,
        title: str,
        minimum: int,
        maximum: int,
        window_modal: bool = True,
        frameless: bool = False,
        size: Optional[Tuple[int, int]] = None,
        cancel_button_text: Optional[str] = None,
        auto_close: bool = False,
        auto_reset: bool = False,
    ) -> QtWidgets.QProgressDialog:
        dialog = QtWidgets.QProgressDialog(label, cancel_button_text, minimum, maximum, self.parent)
        dialog.setWindowTitle(title)
        if window_modal:
            dialog.setWindowModality(Qt.WindowModal)
        if frameless:
            dialog.setWindowFlags(Qt.FramelessWindowHint | Qt.Dialog)
        if size is not None:
            # Use a minimum size instead of a hard fixed size so controls are not clipped
            # under different DPI/font scale settings.
            dialog.setMinimumSize(size[0], size[1])
            dialog.resize(size[0], size[1])
        dialog.setAutoClose(auto_close)
        dialog.setAutoReset(auto_reset)
        # Apply stylesheet BEFORE showing to ensure proper rendering
        if self.stylesheet:
            dialog.setStyleSheet(self.stylesheet)
            # Force stylesheet application to child widgets
            for child in dialog.findChildren(QtWidgets.QWidget):
                if isinstance(child, (QtWidgets.QLabel, QtWidgets.QProgressBar, QtWidgets.QPushButton)):
                    child.setStyleSheet(self.stylesheet)
        dialog.show()
        dialog.setValue(minimum)
        # Process events immediately to ensure dialog is rendered
        QtWidgets.QApplication.processEvents()

        # Ensure cancel button remains fully visible when present.
        cancel_button = dialog.findChild(QtWidgets.QPushButton)
        if cancel_button is not None:
            cancel_button.setMinimumWidth(96)
            cancel_button.setMinimumHeight(30)

        try:
            content_size = dialog.sizeHint()
            target_width = max(dialog.minimumWidth(), content_size.width() + 20)
            target_height = max(dialog.minimumHeight(), content_size.height() + 20)
            dialog.resize(target_width, target_height)
        except Exception:
            pass

        QtWidgets.QApplication.processEvents()
        self.dialog = dialog
        return dialog

    def show_indeterminate(
        self,
        label: str = "Work in progress...",
        title: str = "Progress",
        window_modal: bool = True,
        frameless: bool = False,
        size: Optional[Tuple[int, int]] = None,
        cancel_button_text: Optional[str] = None,
    ) -> QtWidgets.QProgressDialog:
        return self.show(
            label=label,
            title=title,
            minimum=0,
            maximum=0,
            window_modal=window_modal,
            frameless=frameless,
            size=size,
            cancel_button_text=cancel_button_text,
            auto_close=False,
            auto_reset=False,
        )

    def set_label(self, label: str, process_events: bool = False) -> None:
        if self.dialog is None:
            return
        self.dialog.setLabelText(label)
        if process_events:
            QtWidgets.QApplication.processEvents()

    def set_value(self, value: int, process_events: bool = False) -> None:
        if self.dialog is None:
            return
        self.dialog.setValue(value)
        if process_events:
            QtWidgets.QApplication.processEvents()

    def close(self, process_events: bool = True) -> None:
        """Close the progress dialog. Default process_events=True to ensure cleanup."""
        if self.dialog is None:
            return
        self.dialog.hide()
        self.dialog.accept()
        self.dialog.close()
        self.dialog = None
        if process_events:
            QtWidgets.QApplication.processEvents()

    def cancel(self, process_events: bool = True) -> None:
        """Cancel and close the progress dialog. Default process_events=True to ensure cleanup."""
        if self.dialog is None:
            return
        self.dialog.cancel()
        self.dialog = None
        if process_events:
            QtWidgets.QApplication.processEvents()
