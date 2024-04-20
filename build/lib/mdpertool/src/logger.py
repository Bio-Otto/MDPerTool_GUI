import sys
import logging
import logging.config


def Logger(file_name):
    """
    Creates and configures a logger with both file and console handlers.

    :param file_name: The name of the log file.
    :return: The configured logger object.
    """

    # Define the log format
    log_format = '%(asctime)s %(module)s,line: %(lineno)d %(levelname)8s | %(message)s'
    date_format = '%Y/%m/%d %H:%M:%S'

    # Configure the formatter
    formatter = logging.Formatter(fmt=log_format, datefmt=date_format)

    # Configure the file handler
    file_handler = logging.FileHandler(filename=file_name, mode='a')
    file_handler.setFormatter(formatter)

    # Configure the console handler
    console_handler = logging.StreamHandler(stream=sys.stdout)
    console_handler.setFormatter(formatter)

    # Create and configure the logger
    logger = logging.getLogger()
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)
    logger.setLevel(logging.INFO)  # Set the logging level to INFO

    # Add a custom logging level
    logging.addLevelName(35, "DECOMPOSE")
    logging.addLevelName(36, "PREPARATION")
    logging.addLevelName(37, "REFERENCE_MD")
    logging.addLevelName(38, "PERTURBATION_MD")

    # Log a message indicating that the logger object was created successfully
    logger.info("Logger object created successfully.")

    return logger
