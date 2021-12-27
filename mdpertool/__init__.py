import sys
from ._version import __version__, __author__, __credits__, __long_description__, __description__, __url__, \
    __author_email__

# from ui_main import *
if sys.version_info[0] < 3 and 7 > sys.version_info[1] > 9:
    raise Exception('Python 3.7 or greater is required to use this program.')

# from .gui import *
# from .gui.ui_styles import *
# from .src import *
# from .src.omm_runner import *
# from .src.ui_functions import *
# from .src.builder import *
# from .src.mplwidget import *
# from .src.pyside_dynamic import *
# from .src.app_functions import *
# from .ui_main import *
