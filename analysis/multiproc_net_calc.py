from PySide2.QtGui import *
from PySide2.QtWidgets import *
from PySide2.QtCore import *

import time
import traceback, sys


class WorkerSignals(QObject):

    finished = Signal()
    error = Signal(tuple)
    result = Signal(list)
    progress_on_net_calc = Signal(list)
    work_started = Signal()


class Calc_Net_Worker(QRunnable):
    '''
    Worker thread
    Inherits from QRunnable to handler worker thread setup, signals and wrap-up.
    :param callback: The function callback to run on this worker thread. Supplied args and
                     kwargs will be passed through to the runner.
    :type callback: function
    :param args: Arguments to pass to the callback function
    :param kwargs: Keywords to pass to the callback function

    '''

    def __init__(self, fn, *args, **kwargs):
        super(Calc_Net_Worker, self).__init__()

        # Store constructor arguments (re-used for processing)
        self.fn = fn
        self.args = args
        self.kwargs = kwargs
        self.signals = WorkerSignals()

        # Add the callback to our kwargs
        self.kwargs['progress_callback'] = self.signals.progress_on_net_calc

    @Slot()
    def run(self):
        '''
        Initialise the runner function with passed args, kwargs.
        '''
        # Retrieve args/kwargs here; and fire processing using them
        try:
            self.signals.work_started.emit()
            net, log = self.fn(*self.args, **self.kwargs)
        except:
            traceback.print_exc()
            exctype, value = sys.exc_info()[:2]
            self.signals.error.emit((exctype, value, traceback.format_exc()))

        else:
            self.signals.result.emit([net, log])  # Return the result of the processing

        finally:
            self.signals.finished.emit()  # Done