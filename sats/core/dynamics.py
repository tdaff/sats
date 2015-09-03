"""
dynamics.py

Generic dynamics routines.
"""

import os
import sys
import random

import quippy.system
from quippy import Atoms
from quippy.dynamicalsystem import DynamicalSystem
from quippy.io import AtomsWriter
from quippy.potential import Potential, Minim

from sats.core.io import castep_write
from sats.ui.log import info, debug


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

def relax_structure(system, potential, relax_positions=True, relax_cell=True):
    """
    Run a geometry optimisation on the structure to find the energy minimum.

    Parameters
    ----------
    system : ase.Atoms
        A system of atoms to run the minimisation on. The structure is
        altered in-place.
    potential : Potential or str
        A quippy Potential object with the desired potential, or a
        potential_str to initialise a new potential.

    Returns
    -------
    minimised_structure : Atoms
        The geometry optimised structure.
    """
    info("Inside minimiser.")

    qsystem = Atoms(system)

    if not isinstance(potential, Potential):
        potential = Potential(potential)
    qsystem.set_calculator(potential)

    minimiser = Minim(qsystem, relax_positions=relax_positions,
                      relax_cell=relax_cell)

    with Capturing(debug_on_exit=True):
        minimiser.run()

    system.set_cell(qsystem.cell)
    system.set_positions(qsystem.positions)
    system.energy = qsystem.get_potential_energy()

    info("Minimiser done.")

    return system


def molecular_dynamics(system, potential, temperature, total_steps=1100000,
                       timestep=1.0, connect_interval=200, write_interval=20000,
                       equilibration_steps=100000, out_of_plane=None,
                       random_seed=None):
    """
    Run very simple molecular dynamics to generate some configurations. Writes
    configurations out as xyz and CASTEP files.
    """

    info("Inside MD.")
    if random_seed is None:
        random_seed = random.SystemRandom().randint(0, 2**63)
    quippy.system.system_set_random_seeds(random_seed)
    info("Quippy Random Seed {0}.".format(random_seed))
    system = Atoms(system)

    # Can take Potential objects, or just use a string
    if not isinstance(potential, Potential):
        potential = Potential(potential)
    system.set_calculator(potential)

    dynamical_system = DynamicalSystem(system)
    with Capturing(debug_on_exit=True):
        dynamical_system.rescale_velo(temperature)

    if out_of_plane is not None:
        # Stop things moving vertically in the cell
        dynamical_system.atoms.velo[3, :] = 0

    base_dir = os.getcwd()
    run_path = '{0}_{1:g}/'.format(system.info['name'], temperature)
    info("Putting files in {0}.".format(run_path))
    os.mkdir(run_path)
    os.chdir(run_path)

    trajectory = 'traj_{0}_{1:g}.xyz'.format(system.info['name'], temperature)
    out = AtomsWriter(trajectory)

    dynamical_system.atoms.set_cutoff(potential.cutoff() + 2.0)
    dynamical_system.atoms.calc_connect()
    potential.calc(dynamical_system.atoms, force=True, energy=True, virial=True)

    structure_count = 0

    # Basic NVE molecular dynamics
    for step_number in range(1, total_steps + 1):
        dynamical_system.advance_verlet1(
            timestep, virial=dynamical_system.atoms.virial)
        potential.calc(dynamical_system.atoms, force=True, energy=True,
                       virial=True)
        dynamical_system.advance_verlet2(
            timestep, f=dynamical_system.atoms.force,
            virial=dynamical_system.atoms.virial)

        # Maintenance of the system
        if not step_number % connect_interval:
            debug("Connect at step {0}".format(step_number))
            dynamical_system.atoms.calc_connect()
            if step_number < equilibration_steps:
                with Capturing(debug_on_exit=True):
                    dynamical_system.rescale_velo(temperature)

        if not step_number % write_interval:
            debug("Write at step {0}".format(step_number))
            # Print goes to captured stdout
            with Capturing(debug_on_exit=True):
                dynamical_system.print_status(
                    epot=dynamical_system.atoms.energy)
                dynamical_system.rescale_velo(temperature)

            if step_number > equilibration_steps:
                out.write(dynamical_system.atoms)
                sp_path = '{0:03d}'.format(structure_count)
                write_filename = '{0}_{1:g}.{2:03d}'.format(
                    system.info['name'], temperature, structure_count)
                os.mkdir(sp_path)
                os.chdir(sp_path)
                castep_write(dynamical_system.atoms,
                             filename=write_filename)
                info("Wrote a configuration {0}.".format(write_filename))
                os.chdir('..')
                structure_count += 1

    out.close()
    os.chdir(base_dir)

    info("MD Done.")
