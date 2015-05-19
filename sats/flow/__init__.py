"""
sats.flow

General framework for creating and running structures.
"""

from sats.ui import options


def runner():
    """
    Main runner function for generating configurations. Takes no arguments as
    it pulls everything from the options packages.

    """
    structure_types = options.get('structure', 'types')
    for structure in structure_types:
        run_structure(structure)


def run_structure(structure_type):
    """
    Run a single structure type to generate multiple configurations.
    :param structure_type:
    :return:
    """


#!/usr/bin/env python

"""
sats

Frontend to the sats package.

"""

import sys

from quippy import Atoms
from quippy.potential import Potential, Minim

from sats.core.dynamics import molecular_dynamics
from sats.core.io import castep_write
from sats.configurations.sheet import disc
from sats.configurations.sheet import sheet

fs_lattice_parameter = 2.6 # 3.18
fcc_lattice_parameter = 3.95
hcp_lattice_parameter = 2.79
species = 74

sheet_type = sys.argv[1]
angle = float(sys.argv[2])
temperature = int(sys.argv[4])

fs_potential = Potential('IP FS')
fs_potential.set_calc_args({'E_scale': 0.99519, 'r_scale': 0.99302})

lattice_parameter = fs_lattice_parameter

# 6  -> 3.5
# 14 -> 5.2
# 21 -> 6.4
# 30 -> 7.5

if sheet_type == 'disc':
    radius = float(sys.argv[3])
    structure_2d = Atoms(disc(species, radius=radius,
                              alat=fs_lattice_parameter, gamma=angle),
                         pbc=(False, False, False))
    structure_2d.name = "W_{0}_{1:.0f}_{2}".format(sheet_type, angle, radius)
    relax_cell = False

else:
    supercell = int(sys.argv[3])
    structure_2d = Atoms(sheet(species, fs_lattice_parameter, gamma=angle,
                               supercell=supercell, orthorhombic=True))
    structure_2d.name = "W_{0}_{1:.0f}_{2}x{2}".format(sheet_type, angle,
                                                       supercell)
    relax_cell = True

# structure_2d.translate((0, 0, -structure_2d.get_center_of_mass()[2]))


# castep_write(structure_2d, structure_2d.name, optimise=False, fix_lattice=False)

structure_2d.set_cutoff(fs_potential.cutoff() + 2.0)
structure_2d.set_calculator(fs_potential)

minimiser = Minim(structure_2d, relax_positions=True, relax_cell=relax_cell)
minimiser.run()

castep_write(structure_2d, structure_2d.name, optimise=False,
             fix_lattice=True)

# From 10.1007/BF01184339, thermal expansion of tungsten is small, but for
# different temperatures we can expand a bit:
# 3000 K V/V0 = 1.11;
# 5000 K V/V0 = 1.29
#
if temperature == 1000:
    structure_2d.set_lattice(structure_2d.lattice*1.003, scale_positions=True)
elif temperature == 3000:
    structure_2d.set_lattice(structure_2d.lattice*(1.11**(1.0/3)), scale_positions=True)
elif temperature == 5000:
    structure_2d.set_lattice(structure_2d.lattice*(1.29**(1.0/3)), scale_positions=True)

molecular_dynamics(structure_2d, fs_potential, temperature, total_steps=110000,
                   write_interval=2000, equilibration_steps=10000)
