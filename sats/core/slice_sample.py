"""
slice_sample.py

Generic MCMC 'slice sample' routines.
"""

import os
import math
from random import SystemRandom, seed, uniform, randint

import quippy
import quippy.system
from quippy import Atoms
from quippy.io import AtomsWriter

from quippy.potential import Potential

from sats.core.io import castep_write
from sats.ui.log import info, debug


class ParamWriter(object):
    """
    Write out configurations to a trajectory and castep files to
    individual subdirectories.

    """

    def __init__(self, atoms, name, potential, calc_energy=True):
        """
        Initialise the writer with the templating atoms and
        potential.
        """
        self.atoms = atoms
        self.name = name
        self.potential = potential
        self.calc_energy = calc_energy
        self.counter = 0

        # TODO: check naming
        self.base_dir = os.getcwd()
        run_path = '{0}'.format(atoms.info['name'])
        info("Putting files in {0}.".format(run_path))
        os.mkdir(run_path)
        os.chdir(run_path)

        trajectory = 'traj_{0}.xyz'.format(atoms.info['name'])
        self.out = AtomsWriter(trajectory)

    def write_config(self, params):
        cell_params = params[:6]
        # Make a list of 3 item lists
        atom_params = [self.atoms[0].position] + zip(*[iter(params[6:])] * 3)

        new_atoms = self.atoms.copy()
        new_atoms.set_cell([[cell_params[0], 0.0, 0.0],
                            [cell_params[1], cell_params[2], 0.0],
                            [cell_params[3], cell_params[4], cell_params[5]]])

        new_atoms.set_positions(atom_params)
        info("Cell volume: {0}".format(new_atoms.get_volume()))
        if self.calc_energy:
            self.potential.calc(new_atoms, energy=True, forces=True,
                                virial=True)
        # Add to the xyz
        self.out.write(new_atoms)
        sp_path = '{0:03d}'.format(self.counter)
        write_filename = '{0}.{1:03d}'.format(self.atoms.info['name'],
                                              self.counter)
        os.mkdir(sp_path)
        os.chdir(sp_path)
        castep_write(new_atoms, filename=write_filename, kpoint_spacing=0.015)
        info("Wrote a configuration {0}.".format(write_filename))
        os.chdir('..')
        self.counter += 1

    def close(self):
        """
        Clean up.
        """
        self.out.close()
        os.chdir(self.base_dir)


def create_density_function(atoms, potential, pressure, temperature, e0=None):
    """
    Create a function that calculates the density under the given conditions
    that requires only a set of parameters. Everything else is contained
    within the closure.

    Parameters
    ----------
    atoms : ase.Atoms
        Template structure used as the basis of the density function
        calculation.
    potential : Potential or str
        A quippy potential or a potential_str to initialise a potential
        used for calculating the density function.
    pressure : float
        Pressure in eV/A^3 of the state point to calcualte the density function.
    temperature : float
        Temperature in Kelvin of the state point.
    e0 : float
        Energy of the minimised structure. This will be calculated from
        `atoms` if not given.

    """

    if e0 is None:
        potential.calc(atoms, energy=True)
        e0 = atoms.energy

    info("Zero energy: {0}.".format(e0))

    def density_function(params):
        cell_params = params[:6]
        # Make a list of 3 item lists
        atom_params = [atoms[0].position] + zip(*[iter(params[6:])] * 3)

        new_atoms = atoms.copy()
        new_atoms.set_cell([[cell_params[0], 0.0, 0.0],
                            [cell_params[1], cell_params[2], 0.0],
                            [cell_params[3], cell_params[4], cell_params[5]]])

        new_atoms.set_positions(atom_params)

        potential.calc(new_atoms, energy=True)

        # Boltzmann in eV/Kelvin
        value = math.exp(
            -(new_atoms.energy - e0 + pressure * new_atoms.get_volume())
            / (temperature * 0.00008617343))
        return value

    return density_function


def increment_params(density_function, params, idx, delta, max_steps, df_0=None):
    """
    Use slice sampling to generate the next set of parameters changeing the
    value at index.

    Parameters
    ----------
    density_function : function
        A function that takes the parameters as an arguments and returns
        a value of df for those parameters.
    params : list of float
        Initial parameters for the density function.
    idx : integer
        Index of the parameter to change.
    delta : float
        Increment by which to change the value.
    max_steps : int
        Maximum number of steps to increment the value during slicing.
    df_0 : float
        Value of the density function using the input params

    Returns
    -------
    (params, df) : (list of float, float)
        New parameters and the value of the density function at those
        parameters.

    """

    if df_0 is None:
        df_0 = density_function(params)

    params_new = params[:]

    # Level of the slice
    df_slice = uniform(0, 1)*df_0
    slice_l = params[idx] - uniform(0, 1)*delta
    slice_r = slice_l + delta

    slice_increments_left = randint(0, max_steps)

    # widen the slice to the left until it is no longer above df_slice
    # or run out of steps
    for increment in range(slice_increments_left):
        params_new[idx] = slice_l
        df_new = density_function(params_new)
        if df_new > df_slice:
            slice_l -= delta
        else:
            break

    # use the rest of the steps on the right
    for increment in range(max_steps-slice_increments_left):
        params_new[idx] = slice_r
        df_new = density_function(params_new)
        if df_new > df_slice:
            slice_r += delta
        else:
            break

    # Find somewhere in the middle that's good
    while True:
        # Get a point in the middle
        slice_mid = slice_l + uniform(0, 1)*(slice_r-slice_l)
        params_new[idx] = slice_mid
        df_new = density_function(params_new)
        if df_new > df_slice:
            return params_new, df_new
        elif slice_r - slice_r < delta/max_steps:
            # Did it get stuck somewhere and can't get out?
            return params_new, df_new
        elif slice_mid < params[idx]:
            slice_l = params_new[idx]
        else:
            slice_r = params_new[idx]


def slice_sample(bulk, potential, potential_filename, temperature, pressure,
                 lattice_delta, atom_delta, m_max, e0=None, init_d=0,
                 num_configs=10, write_interval=-1, random_seed=None):

    info("Inside Slice Sample.")

    # Randomise the random seed
    if random_seed is None:
        random_seed = SystemRandom().randint(0, 2**63)
    quippy.system.system_set_random_seeds(random_seed)
    seed(random_seed)
    info("Quippy Random Seed {0}.".format(random_seed))
    info("Python Random Seed {0}.".format(random_seed))

    if not isinstance(potential, Potential):
        if potential_filename:
            potential = Potential(potential, param_filename=potential_filename)
        else:
            potential = Potential(potential)

    bulk = Atoms(bulk)
    # pressure in GPa -> eV/A^3
    pressure = pressure/quippy.GPA
    density_function = create_density_function(bulk, potential, pressure,
                                               temperature, e0)

    # Convert to triangle lattice representation (lower triangle in cell)
    scaled_positions = bulk.get_scaled_positions()
    lattice_params = quippy.get_lattice_params(bulk.lattice)
    new_cell = quippy.make_lattice(*lattice_params).T
    bulk.set_cell(new_cell)
    bulk.set_scaled_positions(scaled_positions)
    params = [bulk.cell[0][0], bulk.cell[1][0], bulk.cell[1][1],
              bulk.cell[2][0], bulk.cell[2][1], bulk.cell[2][2]]
    info("Re-orineted to cell: {0}.".format(new_cell.tolist()))

    for atom in bulk.positions.tolist()[1:]:
        params.extend(atom)

    # value for the first iteration
    df_value = density_function(params)

    # Floating count so that division by write_interval make it integer for
    # written configurations
    count = 0.0
    idx = init_d
    output = ParamWriter(bulk, bulk.info['name'], potential)

    # Only write once everything has changed
    if write_interval < 1:
        write_interval = len(params)
    info("Writing configurations after {0} steps.".format(write_interval))

    # Loop just iterates to the next value of x
    while count/write_interval < num_configs:
        # Determine the delta value outside the incrementer, so we can make it
        # do less work
        if idx < 6:
            delta = lattice_delta
        else:
            delta = atom_delta
        # Here's where the magic happens
        params, df_value = increment_params(density_function, params, idx,
                                            delta, m_max, df_0=df_value)

        debug("SLICE_SAMPLE: {0:g}".format(count/write_interval))
        debug("Params: {0}.".format(", ".join("{0}".format(x) for x in params)))

        if not count%write_interval:
            output.write_config(params)

        count += 1
        idx += 1
        if idx >= len(params):
            idx = 0

    output.close()
    info("Slice Sample Done.")
