"""
slice_sample.py

Generic MCMC 'slice sample' routines.
"""

import os
import math
from random import uniform

import quippy
from numpy import dot
from numpy.linalg import inv
from quippy import Atoms
from quippy.dynamicalsystem import DynamicalSystem
from quippy.io import AtomsWriter

from quippy.potential import Potential, Minim

from sats.core.io import castep_write
from sats.ui.log import info, debug
from sats.util import Capturing


def create_density_function(atoms, potential, pressure, temperature, e0=None):
    """
    Create a function that calculates the density under the given conditions
    that requires only a set of parameters. Everything else is contained
    within the closure.

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
        new_atoms.cell = [[cell_params[0], cell_params[1], cell_params[3]],
                          [0.0, cell_params[2], cell_params[4]],
                          [0.0, 0.0, cell_params[5]]]

        new_atoms.set_positions(atom_params)

        potential.calc(new_atoms, energy=True)

        # Boltzmann in eV/Kelvin
        value = math.exp(
            -(new_atoms.energy - e0 + pressure * new_atoms.get_volume())
            / (temperature * 0.00008617343))
        return value

    return density_function


def next_x(density_function, initial_parameters, idx, w_lattice, w_atom, m, f_0=None):
    # make a copy of the list
    x_new = initial_parameters[:]
    x = initial_parameters

    if f_0 is not None:
        f_slice = uniform(0, 1) * f_0
    else:
        f_slice = uniform(0, 1) * density_function(initial_parameters)

    if idx < 6:
        x_new_L = initial_parameters[idx] - (w_lattice * uniform(0, 1))
        x_new_R = x_new_L + w_lattice
    else:
        x_new_L = initial_parameters[idx] - (w_atom * uniform(0, 1))
        x_new_R = x_new_L + w_atom

    J = math.floor(m * uniform(0, 1))
    K = m - 1 - J

    f_x_new_L = 0.0
    f_x_new_R = 0.0
    f_x_new_L_up_to_date = False
    f_x_new_R_up_to_date = False

    while J > 0:
        increase_slice = False

        if not f_x_new_L_up_to_date:
            x[idx] = x_new_L
            f_x_new_L = density_function(x)
            f_x_new_L_up_to_date = True

        if f_x_new_L > f_slice:
            increase_slice = True

        if increase_slice:
            if idx < 7:
                x_new_L = x_new_L - w_lattice
            else:
                x_new_L = x_new_L - w_atom
            f_x_new_L_up_to_date = False
            J -= J
        else:
            break

    while K > 0:
        increase_slice = False

        if not f_x_new_R_up_to_date:
            x[idx] = x_new_R
            f_x_new_R = density_function(x)
            f_x_new_R_up_to_date = True

        if f_x_new_R > f_slice:
            increase_slice = True

        if increase_slice:
            if idx < 7:
                x_new_R = x_new_R + w_lattice
            else:
                x_new_R = x_new_R + w_atom
            f_x_new_R_up_to_date = False
            K -= 1
        else:
            break

    while True:
        x_new[idx] = x_new_L + (uniform(0, 1) * (x_new_R - x_new_L))

        f_x_new = density_function(x_new)

        if f_x_new > f_slice:
            break

        if idx < 7 and x_new_R - x_new_L < w_lattice:
            break
        elif x_new_R - x_new_L < w_atom:
            break
        elif x_new[idx] < initial_parameters[idx]:
            x_new_L = x_new[idx]
        else:
            x_new_R = x_new[idx]

    return x_new, f_x_new

def write_config(atoms, params, name, potential):
    cell_params = params[:6]
    # Make a list of 3 item lists
    atom_params = [atoms[0].position] + zip(*[iter(params[6:])] * 3)

    new_atoms = atoms.copy()
    new_atoms.cell = [[cell_params[0], cell_params[1], cell_params[3]],
                      [0.0, cell_params[2], cell_params[4]],
                      [0.0, 0.0, cell_params[5]]]

    new_atoms.positions = atom_params
    print(new_atoms.get_volume())
    potential.calc(new_atoms, energy=True, forces=True, virial=True)

    new_atoms.write(name)


def slice_sample(bulk, potential, temperature, pressure, lattice_delta,
                 atom_delta, m_max, e0=None, init_d=0, num_configs=10,
                 write_interval=1):

    if not isinstance(potential, Potential):
        potential = Potential(potential)

    bulk = Atoms(bulk)
    # FIXME pressure in GPa? ->
    pressure = pressure/quippy.GPA
    density_function = create_density_function(bulk, potential, pressure,
                                               temperature, e0)

    # Convert to triangle lattice representation
    # TODO: make more transparent
    latt_inv = inv(bulk.lattice)
    lattice_params = quippy.get_lattice_params(bulk.lattice)
    new_lattice = quippy.make_lattice(*lattice_params)
    bulk.lattice = new_lattice
    bulk.positions = [dot(new_lattice, dot(latt_inv, pos)) for pos in
                      bulk.positions]

    new_lattice = new_lattice.tolist()
    params = [new_lattice[0][0], new_lattice[0][1], new_lattice[1][1],
              new_lattice[0][2], new_lattice[1][2], new_lattice[2][2]]

    for atom in bulk.positions.tolist()[1:]:
        params.extend(atom.position)

    df_value = density_function(params)
    print(df_value)

    n = 1
    idx = init_d

    # Loop just iterates to the next value of x
    while n < num_configs*write_interval:
        params, df_value = next_x(density_function, params, idx, lattice_delta, atom_delta,
                   m_max, f_0=df_value)

        info("SLICE_SAMPLE: {0}".format(n))
        if not n%write_interval:
            write_config(bulk, params, 'slice_sample_{0}.xyz'.format(n), potential)

        n += 1
        idx += 1
        if idx >= len(params):
            idx = 0

    info("Slice Sample Done.")
