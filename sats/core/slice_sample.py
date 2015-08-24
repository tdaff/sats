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


def create_density_function(atoms, e0, potential, pressure, temperature):
    """
    Create a function that calculates the density under the given conditions
    that requires only a set of parameters. Everything else is contained
    within the closure.

    """

    def density_function(params):
        cell_params = params[:6]
        # Make a list of 3 item lists
        atom_params = [atoms[0].position] + zip(*[iter(params[6:])] * 3)

        new_atoms = atoms.copy()
        new_atoms.cell = [[cell_params[0], cell_params[1], cell_params[3]],
                          [0.0, cell_params[2], cell_params[4]],
                          [0.0, 0.0, cell_params[5]]]

        new_atoms.positions = atom_params

        potential.calc(new_atoms, energy=True)

        # where did this value come from?
        value = math.exp(
            -(new_atoms.energy - e0 + pressure * new_atoms.get_volume())
            / (temperature * 0.00008617343))
        return value

    return density_function


def next_x(density_function, x_0, d, w_1, w_2, m, f_0=None):
    # make a copy of the list
    x_new = x_0[:] + [0.0]
    x = x_0

    if f_0 is not None:
        f_slice = uniform(0, 1) * f_0
    else:
        f_slice = uniform(0, 1) * density_function()

    if d < 6:
        x_new_L = x_0[d] - (w_1 * uniform(0, 1))
        x_new_R = x_new_L + w_1
    else:
        x_new_L = x_0[d] - (w_2 * uniform(0, 1))
        x_new_R = x_new_L + w_2

    J = math.floor(m * uniform(0, 1))
    K = m - 1 - J

    f_x_new_L = 0.0
    f_x_new_R = 0.0
    f_x_new_L_up_to_date = False
    f_x_new_R_up_to_date = False

    while J > 0:
        increase_slice = False

        if not f_x_new_L_up_to_date:
            x[d] = x_new_L
            f_x_new_L = density_function(x)
            f_x_new_L_up_to_date = True

        if f_x_new_L > f_slice:
            increase_slice = True

        if increase_slice:
            if d < 7:
                x_new_L = x_new_L - w_1
            else:
                x_new_L = x_new_L - w_2
            f_x_new_L_up_to_date = False
            J -= J
        else:
            break

    while K > 0:
        increase_slice = False

        if not f_x_new_R_up_to_date:
            x[d] = x_new_R
            f_x_new_R = density_function(x)
            f_x_new_R_up_to_date = True

        if f_x_new_R > f_slice:
            increase_slice = True

        if increase_slice:
            if d < 7:
                x_new_R = x_new_R + w_1
            else:
                x_new_R = x_new_R + w_2
            f_x_new_R_up_to_date = False
            K -= 1
        else:
            break

    while True:
        x_new[d] = x_new_L + (uniform(0, 1) * (x_new_R - x_new_L))

        f_x_new = density_function(x_new[:-1])

        if f_x_new > f_slice:
            x_new[-1] = f_x_new
            break

        if x_new[d] < x_0[d]:
            x_new_L = x_new[d]
        else:
            x_new_R = x_new[d]

    return x_new

def write_config(atoms, params, name, potential):
    cell_params = params[:6]
    # Make a list of 3 item lists
    atom_params = [atoms[0].position] + zip(*[iter(params[6:])] * 3)

    new_atoms = atoms.copy()
    new_atoms.cell = [[cell_params[0], cell_params[1], cell_params[3]],
                      [0.0, cell_params[2], cell_params[4]],
                      [0.0, 0.0, cell_params[5]]]

    new_atoms.positions = atom_params

    potential.calc(new_atoms, energy=True, forces=True, virial=True)

    new_atoms.write(name)



def slice_sample(bulk, potential, pressure, temperature, e0, lattice_delta, atom_delta, m_max, init_d=0):
    if not isinstance(potential, Potential):
        potential = Potential(potential)

    bulk = Atoms(bulk)
    # FIXME pressure in GPa? ->
    pressure = pressure/quippy.GPA
    density_function = create_density_function(bulk, e0, potential, pressure,
                                               temperature)

    latt_inv = inv(bulk.lattice)
    lattice_params = quippy.get_lattice_params(bulk.lattice)
    new_lattice = quippy.make_lattice(*lattice_params)
    bulk.lattice = new_lattice
    bulk.positions = [dot(new_lattice, dot(latt_inv, x)) for x in
                      bulk.positions]

    new_lattice = new_lattice.tolist()
    x = [new_lattice[0][0], new_lattice[0][1], new_lattice[1][1],
         new_lattice[0][2], new_lattice[1][2], new_lattice[2][2]]

    for atom in bulk[1:]:
        x.extend(atom.position)

    x.append(density_function(x[:-1]))

    n = 1
    d = init_d
    n_configs = 100

    while n < n_configs:
        x_dash = x
        x = next_x(density_function, x_dash[:-1], d, lattice_delta, atom_delta,
                   m_max, f_0=x_dash[-1])

        info("SLICE_SAMPLE: {0}".format(n))
        write_config(bulk, x, 'slice_sample_{0}.xyz'.format(n), potential)

        n += 1
        d += 1
        # includes energy as last item
        if d >= len(x) - 1:
            d = 0

    info("Slice Sample Done.")
