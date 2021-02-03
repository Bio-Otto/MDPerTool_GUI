# import time and OS modules to use to build file folder name
import time
import os
from datetime import datetime
import sys
import logging


def write_folder(directory):
    # Build string for directory to hold files
    # Output Configuration
    folder_name = datetime.now().strftime("%Y-%m-%d_%I-%M-%S_%p")
    created_file_path = os.path.join(directory, folder_name)

    # IF no such folder exists, create one automatically
    if not os.path.exists(created_file_path):
        os.mkdir(created_file_path)

    return created_file_path, folder_name


def lof_file_settings(file_name):
    logging_config = dict(
        version=1,
        formatters={
            'verbose': {
                'format': ("[%(asctime)s] %(levelname)s "
                           "[%(name)s:%(lineno)s] %(message)s"),
                'datefmt': "%d/%b/%Y %H:%M:%S",
            },
            'simple': {
                'format': '%(levelname)s %(message)s',
            },
        },
        handlers={
            'api-logger': {'class': 'logging.handlers.RotatingFileHandler',
                           'formatter': 'verbose',
                           'level': logging.DEBUG,
                           'filename': '%s/api.log' % file_name,
                           'maxBytes': 52428800,
                           'backupCount': 7},
            'batch-process-logger': {'class': 'logging.handlers.RotatingFileHandler',
                                     'formatter': 'verbose',
                                     'level': logging.DEBUG,
                                     'filename': '%s/batch.log' % file_name,
                                     'maxBytes': 52428800,
                                     'backupCount': 7},
            'console': {
                'class': 'logging.StreamHandler',
                'level': 'DEBUG',
                'formatter': 'simple',
                'stream': sys.stdout,
            },
        },
        loggers={
            'api_logger': {
                'handlers': ['api-logger', 'console'],
                'level': logging.DEBUG
            },
            'batch_process_logger': {
                'handlers': ['batch-process-logger', 'console'],
                'level': logging.DEBUG
            }
        }
    )

    return logging_config
