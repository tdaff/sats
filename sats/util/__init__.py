"""
util

General utility functions with no home.

"""

import os
import sys

from sats.ui.log import debug


# CAPTURE quippy output
class Capturing(list):
    """Capture normal stdout and put in a pipe."""
    def __init__(self, debug_on_exit=False, *args, **kwargs):
        self.debug_on_exit = debug_on_exit
        super(Capturing, self).__init__(*args, **kwargs)

    def __enter__(self):
        """Capture normal stdout and put in a pipe."""
        # Use __stdout__ as this works for IPython too
        self.stdout_fileno = sys.__stdout__.fileno()
        self.stdout_save = os.dup(self.stdout_fileno)
        self.pipe_in, self.pipe_out = os.pipe()
        os.dup2(self.pipe_out, self.stdout_fileno)
        os.close(self.pipe_out)
        return self

    def update(self):
        """Read any output so far and update the object."""
        lines = os.read(self.pipe_in, 99999).splitlines()
        self.extend(lines)
        return lines

    def __exit__(self, *args):
        os.close(self.stdout_fileno)
        self.extend(os.read(self.pipe_in, 99999).splitlines())
        os.close(self.pipe_in)
        os.dup2(self.stdout_save, self.stdout_fileno)
        os.close(self.stdout_save)
        if self.debug_on_exit:
            for line in self:
                debug(line)
