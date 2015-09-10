"""
util

General utility functions with no home.

"""

import os
import sys

from numpy.linalg import inv

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


def kpoint_spacing_to_mesh(structure, density):
    """
    Calculate the kpoint mesh that is equivalent to the given density
    in reciprocal Angstrom.

    Parameters
    ----------
    structure : ase.Atoms
        A structure that can have get_reciprocal_cell called on it.
    density : float
        K-Point density in $A^{-1}$.

    Returns
    -------
    kpoint_mesh : [int, int, int]

    """
    # No factor of 2pi in ase, otherwise need to divide through in the mesh
    try:
        r_cell = structure.get_reciprocal_cell()
    except NameError:
        r_cell = inv(structure.cell).transpose()

    r_x = sum(x**2 for x in r_cell[0])**0.5
    r_y = sum(x**2 for x in r_cell[1])**0.5
    r_z = sum(x**2 for x in r_cell[2])**0.5

    kpoint_mesh = [
        int(r_x/(density)) + 1,
        int(r_y/(density)) + 1,
        int(r_z/(density)) + 1]

    debug("Kpoints: {}".format(kpoint_mesh))

    return kpoint_mesh
