import sys
from ._version import __version__, __author__, __credits__, __long_description__, __description__, __url__, __author_email__
#from ui_main import *
if sys.version_info[0] < 3 and 7 > sys.version_info[1] > 8:
    raise Exception('Python 3.7 or greater is required to use this program.')