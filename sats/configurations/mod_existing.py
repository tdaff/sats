"""
mod_existing.py

Helpers modify structures from files.

"""

import re
from random import uniform

import ase.io
from ase import Atom
from numpy import array

from sats.core.elements import Element


def modify(filename, modification, **kwargs):
    """
    Helper function to call the modification of structure contained in
    the file. Use filename@idx to choose an specific frame.

    Parameters
    ----------
    filename : str
        Name of the file to open in the form file.xyz[@idx].
    modification
        Modification type to perform, possible choices are 'add_species'.
    kwargs
        Arguments passed to the modifying procedure.

    Returns
    -------
    structure : ase.Atoms
        Structure after modification

    """

    mods = {
        'add_species': add_species
    }

    atoms = ase.io.read(filename)
    if 'name' not in atoms.info:
        atoms.info['name'] = re.sub(r'\W', '_', filename)

    # Store the forces and energy before modification, seems to raise
    # NotImplementedError for missing info
    try:
        atoms.arrays['base_force'] = atoms.arrays.pop('force')
    except NotImplementedError:
        pass

    try:
        atoms.info['base_energy'] = atoms.get_potential_energy()
        atoms.get_calculator().reset()
    except NotImplementedError:
        pass

    new_atoms = mods[modification](atoms, **kwargs)

    return new_atoms


def add_species(atoms, species='H', count=1, min_dist=0.8, maxiter=10,
                fix_original=True):
    """
    Parameters
    ----------
    atoms : ase.Atoms
        Initial structure to be modified. Name will be updated to reflect
        modification.
    species : str or Element
        Chemical type of the atom to add.
    count : int
        Number of atoms to add.
    min_dist : float
        Keep moving the atoms until the added atom is at least this far from
        everything else.
    maxiter : int
        Number of times to attempt to fit within min_dist. Otherwise just use
        the final configuration.
    fix_original : bool
        Set the move_mask on the input structure to be fixed and added atoms
        allowed to move.

    Returns
    -------
    structure : ase.Atoms
        Modified structure.

    """

    #TODO: finish up here

    fixed_count = len(atoms)

    for _idx in range(count):
        atoms += Atom(symbol=Element(species), position=(0, 0, 0))

        for _trial in range(maxiter):
            new_position = [
                uniform(0, atoms.cell[0][0]) + uniform(0, atoms.cell[1][0]) + uniform(0, atoms.cell[2][0]),
                uniform(0, atoms.cell[0][1]) + uniform(0, atoms.cell[1][1]) + uniform(0, atoms.cell[2][1]),
                uniform(0, atoms.cell[0][2]) + uniform(0, atoms.cell[1][2]) + uniform(0, atoms.cell[2][2])]
            atoms[-1].position = new_position

            if min(atoms.get_distances(-1, range(len(atoms)-1))) > min_dist:
                break

    # set move_mask to 1 only on atom that may move
    if fix_original:
        atoms.arrays['move_mask'] = array([0]*(fixed_count) + [1]*count)

    atoms.info['name'] = "{0}_add_{1}{2}".format(
        atoms.info.get('name', 'structure'), count, species)

    return atoms
