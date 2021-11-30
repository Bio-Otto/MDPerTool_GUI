from PySide2.QtWidgets import QMessageBox
from PySide2.QtCore import Qt
from PySide2.QtGui import QIcon


class Message_Boxes:

    def Information_message(self, title, message, style):
        info_msgbox = QMessageBox(QMessageBox.Information, title, title + "\n\n" + message)
        info_msgbox.setIcon(QMessageBox.Information)
        info_msgbox.addButton(QMessageBox.Ok)
        info_msgbox.setStyleSheet(style)
        info_msgbox.setWindowFlags(Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        info_msgbox.exec_()

    def Critical_message(self, title, message, style):
        critic_msgbox = QMessageBox(QMessageBox.Critical, title, title + "\n\n" + message)
        critic_msgbox.setIcon(QMessageBox.Critical)
        critic_msgbox.addButton(QMessageBox.Ok)
        critic_msgbox.setStyleSheet(style)
        critic_msgbox.setWindowFlags(Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        critic_msgbox.exec_()

    def Warning_message(self, title, message, style):
        warning_msgbox = QMessageBox(QMessageBox.Warning, title, title + "\n\n" + message)
        warning_msgbox.setIcon(QMessageBox.Warning)
        warning_msgbox.addButton(QMessageBox.Ok)
        warning_msgbox.setStyleSheet(style)
        warning_msgbox.setWindowFlags(Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        warning_msgbox.exec_()

    def Succesfully_message(self, title, message, style):
        warning_msgbox = QMessageBox(QMessageBox.Warning, title, title + "\n\n" + message)
        warning_msgbox.setIcon(QMessageBox.Information)

        warning_msgbox.setWindowModality(Qt.NonModal)
        warning_msgbox.setWindowOpacity(.8)

        warning_msgbox.addButton(QMessageBox.Ok)
        warning_msgbox.setStyleSheet(style)
        warning_msgbox.setWindowFlags(Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        warning_msgbox.exec_()

    def Question_message(self, title, message, style):
        question_msgbox = QMessageBox(QMessageBox.Question, title, title + "\n\n" + message)
        question_msgbox.setIcon(QMessageBox.Question)
        question_msgbox.addButton(QMessageBox.Yes)
        question_msgbox.addButton(QMessageBox.No)
        question_msgbox.setDefaultButton(QMessageBox.No)
        # msgbox.setCheckBox(cb)
        # msg.setWindowIcon(QIcon(ICON_PATH))
        question_msgbox.setStyleSheet(style)
        question_msgbox.setWindowFlags(Qt.FramelessWindowHint | Qt.WindowStaysOnTopHint)
        close_answer = question_msgbox.exec_()

        return close_answer
