#!/usr/bin/env python

"""
Run sats and generate some configurations.

"""

import sys

import quippy
from quippy.potential import Potential, Minim
from quippy.io import AtomsWriter
from quippy.dynamicalsystem import DynamicalSystem, BAROSTAT_HOOVER_LANGEVIN, THERMOSTAT_ALL_PURPOSE
from quippy.structures import bcc, fcc, supercell


fs_lattice_parameter = 3.18
species = 74
lattice = quippy.bcc

# This is a little bit higher as python sometimes truncates 3.99
max_atoms = 512

if lattice.func_name == 'sc':
    lattice_parameter = fs_lattice_parameter*0.5*(2**0.5)
else:
    lattice_parameter = fs_lattice_parameter

# Make a cubic supercell with up to max_atoms in it
bulk = lattice(lattice_parameter, species)
n_supercell = int((float(max_atoms)/bulk.n)**(1.0/3.0))
bulk = supercell(bulk, n_supercell, n_supercell, n_supercell)

fs_potential = Potential('IP FS')
fs_potential.set_calc_args({'E_scale': 0.99519, 'r_scale': 0.99302})

bulk.set_cutoff(fs_potential.cutoff() + 2.0)
bulk.set_calculator(fs_potential)

minimiser = Minim(bulk, relax_positions=True, relax_cell=True)
minimiser.run()

TEMPERATURE = 5000
# From 10.1007/BF01184339, thermal expansion of tungsten is small, but for
# different temperatures we can expand a bit:
# 3000 K V/V0 = 1.11;
# 5000 K V/V0 = 1.29
#
if TEMPERATURE == 1000:
    bulk.set_lattice(bulk.lattice*1.003, scale_positions=True)
elif TEMPERATURE == 3000:
    bulk.set_lattice(bulk.lattice*(1.11**(1.0/3)), scale_positions=True)
elif TEMPERATURE == 5000:
    bulk.set_lattice(bulk.lattice*(1.29**(1.0/3)), scale_positions=True)

dynamical_system = DynamicalSystem(bulk)
dynamical_system.set_barostat(type=BAROSTAT_HOOVER_LANGEVIN, p_ext=0.0,
                              hydrostatic_strain=True, diagonal_strain=True,
                              finite_strain_formulation=False, tau_epsilon=10.0,
                              w_epsilon=10000.0, w_epsilon_factor=10000.0,
                              t=TEMPERATURE)
dynamical_system.add_thermostat(type=THERMOSTAT_ALL_PURPOSE, t=TEMPERATURE, tau=100.0)
dynamical_system.rescale_velo(TEMPERATURE)

total_steps = 1100000
timestep = 0.1  # fs
connect_interval = 200
write_interval = 20000
equilibration_steps = 100000

trajectory = 'traj_{}_{}.xyz'.format(lattice.func_name, TEMPERATURE)
out = AtomsWriter(trajectory)

dynamical_system.atoms.calc_connect()
fs_potential.calc(dynamical_system.atoms, force=True, energy=True, virial=True)

# Basic NVE molecular dynamics
for step_number in range(1, total_steps+1):
    dynamical_system.advance_verlet1(timestep,
                                     virial=dynamical_system.atoms.virial)
    fs_potential.calc(dynamical_system.atoms, force=True,
                      energy=True, virial=True)
    dynamical_system.advance_verlet2(timestep, f=dynamical_system.atoms.force,
                                     virial=dynamical_system.atoms.virial)

    # Maintenance of the system
    if not step_number % connect_interval:
        dynamical_system.atoms.calc_connect()
        if step_number < equilibration_steps:
            dynamical_system.rescale_velo(TEMPERATURE)

    if not step_number % write_interval:
        dynamical_system.print_status(epot=dynamical_system.atoms.energy)
        dynamical_system.rescale_velo(TEMPERATURE)
        if step_number > equilibration_steps:
            out.write(dynamical_system.atoms)

out.close()
