"""
dynamics.py

Generic dynamics routines.
"""

import os

from quippy.dynamicalsystem import DynamicalSystem
from quippy.io import AtomsWriter

from sats.io import castep_write


def molecular_dynamics(system, potential, temperature, total_steps=1100000,
                       timestep=1.0, connect_interval=200, write_interval=20000,
                       equilibration_steps=100000):

    dynamical_system = DynamicalSystem(system)
    dynamical_system.rescale_velo(temperature)

    base_dir = os.getcwd()
    run_path = '{0}_{1}/'.format(system.name, temperature)
    os.mkdir(run_path)
    os.chdir(run_path)

    trajectory = 'traj_{0}_{1}.xyz'.format(system.name, temperature)
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
                write_filename = '{0}_{1}.{2:03d}'.format(
                    system.name, temperature, structure_count)
                os.mkdir(sp_path)
                os.chdir(sp_path)
                castep_write(dynamical_system.atoms,
                             filename=write_filename)
                os.chdir('..')
                structure_count += 1

    out.close()
    os.chdir(base_dir)
