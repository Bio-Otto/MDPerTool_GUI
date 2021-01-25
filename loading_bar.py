import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QVBoxLayout, QDesktopWidget
from PyQt5.QtCore import Qt, QTimer, QByteArray
from PyQt5.QtGui import QMovie
from PyQt5 import QtCore
from PyQt5.Qt import QSizePolicy


class LoadingScreen(QWidget):

    def __init__(self, filename, title):
        super(LoadingScreen, self).__init__()

        # Load the file into a QMovie
        self.movie = QMovie(filename, QByteArray(), self)

        size = self.movie.scaledSize()
        self.setGeometry(200, 200, size.width(), size.height())
        self.setWindowTitle(title)

        self.movie_screen = QLabel()
        # Make label fit the gif
        self.movie_screen.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.movie_screen.setAlignment(Qt.AlignCenter)

        # Create the layout
        main_layout = QVBoxLayout()
        main_layout.addWidget(self.movie_screen)

        self.setLayout(main_layout)

        # Add the QMovie object to the label
        self.movie.setCacheMode(QMovie.CacheAll)
        self.movie.setSpeed(100)
        self.movie_screen.setMovie(self.movie)
        self.movie.start()

        # timer = QTimer(self)
        # timer.singleShot(3000, self.stop_animation)

    def stop_animation(self):
        self.movie.stop()
        self.close()

    def startAnimation(self):
        self.movie.start()


# if __name__ == "__main__":
#     gif = "C:\\Users\\HIbrahim\\Desktop\\MDPERTOOL_v01\\icons\\gifs\\please_wait.gif"
#     app = QApplication(sys.argv)
#     player = LoadingScreen(gif, "was")
#     player.show()
#     sys.exit(app.exec_())
