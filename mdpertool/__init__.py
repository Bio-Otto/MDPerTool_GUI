"""Top-level package for mdpertool."""

__author__ = """H.Ibrahim Ozdemir"""
__email__ = "halil.ibrahim.oozdemir@gmail.com"
__version__ = '0.0.1'

import os
import PySide2

dirname = os.path.dirname(PySide2.__file__)
plugin_path = os.path.join(dirname, 'plugins', 'platforms')
os.environ['QT_QPA_PLATFORM_PLUGIN_PATH'] = plugin_path

from . import gui
from . import analysis
from mdpertool import *
#from .utils import *
