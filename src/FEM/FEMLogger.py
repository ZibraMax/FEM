#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# -------------------------------------------------------------------------------
#                                                                               -
#  Python dual-logging setup (console and log file),                            -
#  supporting different log levels and colorized output                         -
#                                                                               -
#  Created by Fonic <https://github.com/fonic>                                  -
#  Date: 04/05/20                                                               -
#                                                                               -
#  Based on:                                                                    -
#  https://stackoverflow.com/a/13733863/1976617                                 -
#  https://uran198.github.io/en/python/2016/07/12/colorful-python-logging.html  -
#  https://en.wikipedia.org/wiki/ANSI_escape_code#Colors                        -
#                                                                               -
# -------------------------------------------------------------------------------

# Imports
import os
import sys
import logging
from datetime import datetime


class TimeFilter(logging.Filter):

    def filter(self, record):
        try:
            last = self.last
        except AttributeError:
            last = record.relativeCreated

        delta = datetime.fromtimestamp(
            record.relativeCreated/1000.0) - datetime.fromtimestamp(last/1000.0)

        record.relative = '{0:.4f}'.format(
            delta.seconds + delta.microseconds/1000000.0)

        self.last = record.relativeCreated
        return True


class LogFormatter(logging.Formatter):
    """Creates a Logging Formatter

    Args:
        color (color): Color
    """

    COLOR_CODES = {
        logging.CRITICAL: "\033[1;35m",  # bright/bold magenta
        logging.ERROR:    "\033[1;31m",  # bright/bold red
        logging.WARNING:  "\033[1;33m",  # bright/bold yellow
        logging.INFO:     "\033[0;37m",  # white / light gray
        logging.DEBUG:    "\033[1;30m"  # bright/bold black / dark gray
    }

    RESET_CODE = "\033[0m"

    def __init__(self, color, *args, **kwargs):
        """Creates a Logging Formatter

        Args:
            color (color): Color
        """
        super(LogFormatter, self).__init__(*args, **kwargs)
        self.color = color

    def format(self, record, *args, **kwargs):
        """Formats arecord

        Args:
            record (record): Record to be formatted

        """
        if (self.color == True and record.levelno in self.COLOR_CODES):
            record.color_on = self.COLOR_CODES[record.levelno]
            record.color_off = self.RESET_CODE
        else:
            record.color_on = ""
            record.color_off = ""
        return super(LogFormatter, self).format(record, *args, **kwargs)


class FEMLogger():
    """Creation of a Logger for FEM purposes. Based on Python Logger by Fonic <https://github.com/fonic>   
    """

    def __init__(self, ad: str = '') -> None:
        """Creates a FEMLogger object

        Args:
            ad (str, optional): Aditional name. Defaults to ''.
        """
        self.additional_name = ad

    def setup_logging(self, console_log_output="stdout", console_log_level="warning", console_log_color=True, logfile_file=None, logfile_log_level="debug", logfile_log_color=False, log_line_template="%(color_on)s[%(levelname)-8s] %(message)s%(color_off)s"):
        """Set the logger

        Args:
            console_log_output (str, optional): . Defaults to "stdout".
            console_log_level (str, optional): . Defaults to "warning".
            console_log_color (bool, optional): . Defaults to True.
            logfile_file (optional): . Defaults to None.
            logfile_log_level (str, optional): . Defaults to "debug".
            logfile_log_color (bool, optional): . Defaults to False.
            log_line_template (str, optional): . Defaults to "%(color_on)s[%(levelname)-8s] %(message)s%(color_off)s".

        """
        script_name = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        additional_name = ''
        if self.additional_name:
            additional_name = '_' + self.additional_name
        if not logfile_file:
            logfile_file = script_name + additional_name + ".log"
        # Create logger
        # For simplicity, we use the root logger, i.e. call 'logging.getLogger()'
        # without name argument. This way we can simply use module methods for
        # for logging throughout the script. An alternative would be exporting
        # the logger, i.e. 'global logger; logger = logging.getLogger("<name>")'
        logger = logging.getLogger()
        logging.getLogger('matplotlib').setLevel(logging.WARNING)
        logging.getLogger('PIL').setLevel(logging.WARNING)
        # Set global log level to 'debug' (required for handler levels to work)
        logger.setLevel(logging.DEBUG)

        # Create console handler
        console_log_output = console_log_output.lower()
        if (console_log_output == "stdout"):
            console_log_output = sys.stdout
        elif (console_log_output == "stderr"):
            console_log_output = sys.stderr
        else:
            print("Failed to set console output: invalid output: '%s'" %
                  console_log_output)
            return False
        console_handler = logging.StreamHandler(console_log_output)

        # Set console log level
        try:
            # only accepts uppercase level names
            console_handler.setLevel(console_log_level.upper())
        except Exception as e:
            print("Failed to set console log level: invalid level: '%s'" %
                  console_log_level, e)
            return False

        # Create and set formatter, add console handler to logger
        console_formatter = LogFormatter(
            fmt=log_line_template, color=console_log_color)
        console_handler.setFormatter(console_formatter)
        logger.addHandler(console_handler)

        # Create log file handler
        try:
            logfile_handler = logging.FileHandler(logfile_file)
        except Exception as exception:
            print("Failed to set up log file: %s" % str(exception))
            return False

        # Set log file log level
        try:
            # only accepts uppercase level names
            logfile_handler.setLevel(logfile_log_level.upper())
        except:
            print("Failed to set log file log level: invalid level: '%s'" %
                  logfile_log_level)
            return False

        # Create and set formatter, add log file handler to logger
        logfile_formatter = LogFormatter(
            fmt='[%(asctime)s] (Delta duration: %(relative)ss) ' + log_line_template, color=logfile_log_color)
        logfile_handler.setFormatter(logfile_formatter)
        logger.addHandler(logfile_handler)
        [hndl.addFilter(TimeFilter()) for hndl in logger.handlers]
        self.start_time = datetime.now()
        logging.debug(
            f'============================{script_name}============================')
        logging.debug(
            f'Session started @ {self.start_time.strftime("%d/%m/%Y - %H:%M:%S")}')
        # Success

        return True

    def end_timer(self):
        """Ends the sesion time
        """
        self.end_time = datetime.now()
        logging.debug(
            f'Session ended @ {self.end_time.strftime("%d/%m/%Y - %H:%M:%S")}')
        logging.debug(
            f'Duration: {(self.end_time-self.start_time)}')
        return self.end_time-self.start_time
