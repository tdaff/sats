# -*- coding: utf-8 -*-

"""
log.py

Logging functions that spew output to the logging module

"""

from __future__ import print_function


__all__ = ['debug', 'info', 'warning', 'error', 'critical']

import copy
import logging
import sys
from logging import debug, info, warning, error, critical


class ColouredConsoleHandler(logging.StreamHandler):
    """Makes colourised and wrapped output for the console."""
    def emit(self, record):
        """Colourise and emit a record."""
        # Need to make a actual copy of the record
        # to prevent altering the message for other loggers
        myrecord = copy.copy(record)
        levelno = myrecord.levelno
        if levelno >= 50:  # CRITICAL / FATAL
            front = '\033[30;41m'  # black/red
        elif levelno >= 40:  # ERROR
            front = '\033[30;41m'  # black/red
        elif levelno >= 30:  # WARNING
            front = '\033[30;43m'  # black/yellow
        elif levelno >= 20:  # INFO
            front = '\033[30;42m'  # black/green
        elif levelno >= 10:  # DEBUG
            front = '\033[30;46m'  # black/cyan
        else:  # NOTSET and anything else
            front = '\033[0m'  # normal

        myrecord.levelname = '%s%s\033[0m' % (front, myrecord.levelname)
        logging.StreamHandler.emit(self, myrecord)


def _init_logging(verbosity=0):
    """
    Setup the logging to terminal and sats.log file, with levels as required.

    """
    root_logger = logging.getLogger()

    # Have the logger itself set with the lowest possible level
    root_logger.setLevel(logging.DEBUG)
    # Reset any handlers that might have been set accidentally
    root_logger.handlers = []

    # Always at least INFO in .flog
    file_level = logging.INFO
    log_filename = 'sats.log'

    if verbosity <= -2:
        stdout_level = logging.CRITICAL
    elif verbosity <= -1:
        stdout_level = logging.ERROR
    elif verbosity >= 1:
        stdout_level = logging.DEBUG
        file_level = logging.DEBUG
    else:
        stdout_level = logging.INFO

    # Easier to do simple file configuration then add the stdout
    file_handler = logging.FileHandler(log_filename)
    file_handler.setLevel(file_level)
    formatter = logging.Formatter('[%(asctime)s] %(levelname)s <%(funcName)s> %(message)s',
                                  datefmt='%Y%m%d %H:%M:%S')
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)

    # Make these uniform widths
    logging.addLevelName(10, '--')
    logging.addLevelName(20, '>>')
    logging.addLevelName(30, '**')
    logging.addLevelName(40, '!!')
    logging.addLevelName(50, 'XX')

    # Use nice coloured console output
    console = ColouredConsoleHandler(stream=sys.stdout)
    console.setLevel(stdout_level)
    formatter = logging.Formatter('%(levelname)s %(message)s')
    console.setFormatter(formatter)
    # add the handler to the root logger
    root_logger.addHandler(console)

_init_logging()
