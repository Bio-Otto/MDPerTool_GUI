import sys
from ._version import __version__, __author__, __credits__
#from ui_main import *
if sys.version_info[0] < 3 and sys.version_info[1] < 7:
    raise Exception('Python 3.7 or greater is required to use this program.')