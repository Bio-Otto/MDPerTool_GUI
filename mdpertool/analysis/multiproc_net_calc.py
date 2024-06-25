from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCore import *

import time
import traceback
import sys
from typing import Callable, Any, Tuple, List


class WorkerSignals(QObject):
    """
    Defines the signals available from a running worker thread.
    """
    finished = Signal()
    error = Signal(tuple)
    result = Signal(list)
    progress_on_net_calc = Signal(list)
    work_started = Signal()


class CalcNetWorker(QRunnable):
    """
    Worker thread for running a function in a separate thread.

    Inherits from QRunnable to handle worker thread setup, signals, and wrap-up.

    :param fn: The function callback to run on this worker thread. Supplied args and
               kwargs will be passed through to the runner.
    :type fn: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function
    """
    def __init__(self, fn: Callable, *args: Any, **kwargs: Any):
        super(CalcNetWorker, self).__init__()
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the progress callback to our kwargs
        self.kwargs['progress_callback'] = self.signals.progress_on_net_calc

    @Slot()
    def run(self):
        """
        Initialize the runner function with passed args and kwargs.
        """
        try:
            self.signals.work_started.emit()
            net, log = self.fn(*self.args, **self.kwargs)
            time.sleep(0.2)
        except Exception as e:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))
        else:
            self.signals.result.emit([net, log])
            print(log)
        finally:
            self.signals.finished.emit()

