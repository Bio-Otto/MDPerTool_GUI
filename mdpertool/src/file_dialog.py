from PySide2.QtWidgets import QFileDialog, QDialog
import os


class Dialog():

    def __init__(self, mainform):
        self.__mainform = mainform
        self.__dialog = QFileDialog()
        self.__directory = ''
        self.__filename = ['', '', '']
        self.__filters = []
        self.__default_filter_index = 0
        self.__path = ''

    @property
    def path(self):
        return self.__path

    @property
    def filename(self):
        return self.__filename

    @property
    def directory(self):
        return self.__directory

    @directory.setter
    def directory(self, value):
        self.__directory = value

    @property
    def filters(self):
        return self.__filters

    @filters.setter
    def filters(self, value):
        self.__filters = value

    @property
    def default_filter_index(self):
        return self.__default_filter_index

    @default_filter_index.setter
    def default_filter_index(self, value):
        self.__default_filter_index = value

    def exec(self, load):
        self.__dialog.setNameFilters(self.__filters)
        self.__dialog.selectNameFilter(self.__filters[self.__default_filter_index])
        self.__dialog.setDirectory(self.__directory)
        if load:
            self.__dialog.setLabelText(QFileDialog.Accept, 'Open')
            self.__dialog.setWindowTitle('Open')

        else:
            self.__dialog.setLabelText(QFileDialog.Accept, 'Save')
            self.__dialog.setWindowTitle('Save')

        if self.__dialog.exec() == QDialog.Accepted:
            self.__path = self.__dialog.selectedFiles()[0]
            fn = os.path.split(self.__path)
            ex = os.path.splitext(self.__path)[1]
            self.__filename = [fn[0], fn[1], ex[1:len(ex)]]