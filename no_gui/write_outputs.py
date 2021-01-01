# import time and OS modules to use to build file folder name
import time
import os
from datetime import datetime


def write_folder(directory):
    # Build string for directory to hold files
    # Output Configuration

    folder_name = datetime.now().strftime("%Y-%m-%d_%I-%M-%S_%p")

    created_file_path = os.path.join(directory, folder_name)
    # IF no such folder exists, create one automatically
    if not os.path.exists(created_file_path):
        os.mkdir(created_file_path)

    return created_file_path
