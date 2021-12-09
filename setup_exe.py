"""
# --------------------------------------------------> INFORMATIONS < --------------------------------------------------#
---> DEVELOPED BY: HALİL İBRAHİM ÖZDEMİR
---> PROJECT MADE WITH: PySide2
---> VERSION: 0.1
# --------------------------------------------------> INFORMATIONS < --------------------------------------------------#
"""

import sys
# import os
# import PySide2
from cx_Freeze import setup, Executable

base = None
if sys.platform == "win32":
    base = "Win32GUI"

# dependencies
build_exe_options = {
    "packages": ["os", "sys", "re", "PySide2.QtCore", "PySide2.QtWidgets", "PySide2.QtUiTools", "PySide2.QtQuick",
                 "PySide2.QtQml", "PySide2.QtGui", "matplotlib", "pandas", "pyqtgraph", "OpenGL", "pymol",
                 "networkx", "prody", "pystache", "shiboken2", "pdbfixer", "pyqtgraph", "openmm",
                 "mdtraj", "parmed", "pyvis", "OpenGL"],


    "include_files": ['MAIN_GUI.ui', 'LICENSE', 'splash_screen.ui',
                      'C:\\Users\\HIbrahim\\anaconda3\\envs\\moldyn\\Lib\\xdrlib.py', 'test/'],
    # 'C:\\Users\\HIbrahim\\anaconda3\\envs\\MolDynAnalyze\\Lib\\xdrlib.py',
    "build_exe": "build"
}

target = [
    Executable("ui_main.py",
               base=base,
               target_name="MDPerTool.exe",
               icon="icons/big_icons/style_icon_48x48.png"
               )
]

setup(name="MDPerTool",
      version="0.1",
      description="Perturbation based Allosteric Pathway Finder",
      options={"build_exe": build_exe_options},
      executables=target,
      author="H. Ibrahim Ozdemir",
      )
