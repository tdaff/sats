"""
dynamics.py

Generic dynamics routines.
"""

import os

from quippy import Atoms
from quippy.dynamicalsystem import DynamicalSystem
from quippy.io import AtomsWriter
from quippy.potential import Potential, Minim

from sats.core.io import castep_write


def relax_structure(system, potential, relax_positions=True, relax_cell=True):
    """
    Run a geometry optimisation on the structure to find the energy minimum.

    Parameters
    ----------
    system : ase.Atoms
        A system of atoms to run the minimisation on. The structure is
        altered in-place.
    potential : str
        The potential arg_str to pass to the quippy Potential object.

    Returns
    -------
    minimised_structure : Atoms
        The geometry optimised structure.
    """
    run_potential = Potential(potential)
    system.set_calculator(run_potential)

    minimiser = Minim(system, relax_positions=relax_positions,
                      relax_cell=relax_cell)
    minimiser.run()

    return system


def molecular_dynamics(system, potential_str, temperature, total_steps=1100000,
                       timestep=1.0, connect_interval=200, write_interval=20000,
                       equilibration_steps=100000):
    """
    Run very simple molecular dynamics to generate some configurations. Writes
    configurations out as xyz and CASTEP files.
    """

    system = Atoms(system)

    potential = Potential(potential_str)
    system.set_calculator(potential)

    dynamical_system = DynamicalSystem(system)
    dynamical_system.rescale_velo(temperature)
    print(dynamical_system.atoms.velo.T)
    dynamical_system.atoms.velo[3, :] = 0
    print(dynamical_system.atoms.velo.T)

    base_dir = os.getcwd()
    run_path = '{0}_{1:g}/'.format(system.name, temperature)
    os.mkdir(run_path)
    os.chdir(run_path)

    trajectory = 'traj_{0}_{1:g}.xyz'.format(system.name, temperature)
    out = AtomsWriter(trajectory)

    dynamical_system.atoms.calc_connect()
    potential.calc(dynamical_system.atoms, force=True, energy=True, virial=True)

    structure_count = 0

    # Basic NVE molecular dynamics
    for step_number in range(1, total_steps + 1):
        dynamical_system.advance_verlet1(timestep,
                                         virial=dynamical_system.atoms.virial)
        potential.calc(dynamical_system.atoms, force=True, energy=True,
                       virial=True)
        dynamical_system.advance_verlet2(timestep,
                                         f=dynamical_system.atoms.force,
                                         virial=dynamical_system.atoms.virial)

        # Maintenance of the system
        if not step_number % connect_interval:
            dynamical_system.atoms.calc_connect()
            if step_number < equilibration_steps:
                dynamical_system.rescale_velo(temperature)

        if not step_number % write_interval:
            dynamical_system.print_status(epot=dynamical_system.atoms.energy)
            dynamical_system.rescale_velo(temperature)
            if step_number > equilibration_steps:
                out.write(dynamical_system.atoms)
                sp_path = '{0:03d}'.format(structure_count)
                write_filename = '{0}_{1:g}.{2:03d}'.format(
                    system.name, temperature, structure_count)
                os.mkdir(sp_path)
                os.chdir(sp_path)
                castep_write(dynamical_system.atoms,
                             filename=write_filename)
                os.chdir('..')
                structure_count += 1

    out.close()
    os.chdir(base_dir)
