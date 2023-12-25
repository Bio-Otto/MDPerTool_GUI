import sys
import os
from PySide2.QtWidgets import QApplication, QMainWindow, QTextEdit, QPushButton, QVBoxLayout, QWidget, QFileDialog
from PySide2.QtGui import QColor
import subprocess

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Script Çalıştırma")

        # Metin düzenleyici oluştur
        self.text_edit = QTextEdit()

        # Dosya seç düğmesi oluştur
        self.select_file_button = QPushButton("Dosya Seç")
        self.select_file_button.clicked.connect(self.select_file)

        # Çalıştırma düğmesi oluştur
        self.run_button = QPushButton("Çalıştır")
        self.run_button.clicked.connect(self.run_script)

        # Durdurma düğmesi oluştur
        self.stop_button = QPushButton("Durdur")
        self.stop_button.setEnabled(False)
        self.stop_button.clicked.connect(self.stop_script)

        # Çıktı etiketi oluştur
        self.output_label = QTextEdit()
        self.output_label.setReadOnly(True)

        # Ana düzen oluştur
        layout = QVBoxLayout()
        layout.addWidget(self.text_edit)
        layout.addWidget(self.select_file_button)
        layout.addWidget(self.run_button)
        layout.addWidget(self.stop_button)
        layout.addWidget(self.output_label)

        # Ana widget oluştur ve düzeni ayarla
        widget = QWidget()
        widget.setLayout(layout)
        self.setCentralWidget(widget)

        # Çalışan subprocess nesnesi
        self.subprocess = None

    def select_file(self):
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, "Dosya Seç", "", "Text Files (*.txt)")
        if file_path:
            with open(file_path, 'r') as file:
                script_content = file.read()
                self.text_edit.setPlainText(script_content)

    def run_script(self):
        # Metin düzenleyiciden script içeriğini al
        script_content = self.text_edit.toPlainText()

        # Geçici bir Python dosyası oluştur
        temp_script_file = "temp_script.py"
        with open(temp_script_file, 'w') as file:
            file.write(script_content)

        # Geçici Python dosyasını subprocess ile çalıştır
        try:
            print("yes1")
            self.subprocess = subprocess.Popen([sys.executable, temp_script_file], stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE,
                                               universal_newlines=True)
            print("yes2")

            self.stop_button.setEnabled(True)
            self.run_button.setEnabled(False)

            # Çıktıyı oku ve göster
            while True:
                line = self.subprocess.stdout.readline()
                if not line:
                    break
                if "Error" in line:
                    self.output_label.append(f'<font color="red">{line}</font>')
                else:
                    self.output_label.append(f'<font color="green">{line}</font>')

        except subprocess.CalledProcessError as e:
            error_output = e.output if e.output else str(e)
            self.output_label.append(f'<font color="red">{error_output}</font>')

        # Geçici dosyayı sil
        os.remove(temp_script_file)

        # Durdurma düğmesini devre dışı bırak, diğer düğmeleri etkinleştir
        self.stop_button.setEnabled(False)
        self.run_button.setEnabled(True)

    def stop_script(self):
        # Subprocess'ı sonlandır
        self.subprocess.terminate()

        # Durdurma düğmesini devre dışı bırak, diğer düğmeleri etkinleştir
        self.stop_button.setEnabled(False)
        self.run_button.setEnabled(True)

"""
if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())
"""

import numpy as np
new_equilizer_data = np.column_stack((range(100, 200), np.arange(10, 0, -0.1)))
print(new_equilizer_data)