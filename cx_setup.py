import sys
import os
from cx_Freeze import setup, Executable

base = None
if sys.platform == "win32":
    base = "Win32GUI"

# ADD FILES
files = ['icon.png', 'LICENSE', 'fonts/', 'test/', 'gui/', 'src/', 'analysis/', 'Download/',
         'no_gui/', 'C:\\Users\\law5_\\.conda\\envs\\mdpertool\\Lib\\xdrlib.py']

# TARGET
target = Executable(
    script="ui_main.py",
    base=base,
    icon="icon.png"
)

# PACKAGES
packages = ["OpenGL", "pymol", "parmed", "mdtraj", "openmm"]
#                   "os", "sys", "re", "matplotlib", "pandas", "pyqtgraph", "OpenGL", "pymol", "numpy",
#                  "networkx", "prody", "pystache", "shiboken2", "pdbfixer", "pyqtgraph", "openmm",
#                  "mdtraj", "parmed", "pyvis"

# SETUP CX FREEZE
setup(
    name="MDPerTool",
    version="0.1",
    description="Allosteric Pathway Predicter",
    author="Halil I. Ozdemir",
    options={'build_exe': {'include_files': files, 'packages': packages}},
    executables=[target]

)

# base = None
# if sys.platform == "win32":
#     base = "Win32GUI"
#
# # dependencies
# build_exe_options = {
#     "packages": ["os", "sys", "re", "PySide2.QtCore", "PySide2.QtWidgets", "PySide2.QtUiTools", "PySide2.QtQuick",
#                  "PySide2.QtQml", "PySide2.QtGui", "matplotlib", "pandas", "pyqtgraph", "OpenGL", "pymol", "numpy",
#                  "networkx", "prody", "pystache", "shiboken2", "pdbfixer", "pyqtgraph", "openmm",
#                  "mdtraj", "parmed", "pyvis", "OpenGL"],
#
#     "include_files": ['LICENSE', 'fonts/', 'test/', 'gui/', 'src/', 'analysis/', 'Download/', 'no_gui/'],
#     "excludes": ["tkinter", "PyQt5"],
#     "build_exe": "build"
# }
#
# target = [
#     Executable("ui_main.py",
#                base=base,
#                target_name="MDPerTool.exe",
#                icon="icon.png"
#                )
# ]
#
# setup(name="MDPerTool_GUI",
#       version="0.1",
#       description="Perturbation based Allosteric Pathway Finder",
#       options={"build_exe": build_exe_options},
#       executables=target,
#       author="H. Ibrahim Ozdemir",
#       )
