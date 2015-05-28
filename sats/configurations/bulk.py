"""

Bulk material configurations.

"""

from __future__ import division

from math import ceil

from ase import Atoms
from ase.lattice import bulk as ase_bulk

from sats.bulk.rescale import properties


def bulk(species, lattice, min_atoms=0, max_atoms=None):
    """
    Helper to generate a generic bulk configuration.

    Parameters
    ----------
    species : str or sats.core.elements.Element
        Atomic species used to generate the structures.
    lattice : str
        Bulk lattice type to generate.
    min_atoms : int
        Minimum number of atoms to include in a supercell of the
        bulk.
    max_atoms : int, optional
        Not implemented.

    Returns
    -------
    bulk : ase.Atoms
        Atoms in the desired lattice.
    """

    lattice_constant = properties[species]['lattice_constant'][lattice]

    if lattice == 'a15':
        atoms = a15(lattice_constant, species)
    elif lattice in ['fcc']:
        # Special check for fcc since the cubic and orthorhombic are
        # different.
        atoms = ase_bulk(species, lattice, a=lattice_constant,
                         cubic=True)
    else:
        atoms = ase_bulk(species, lattice, a=lattice_constant,
                         orthorhombic=True)

    # TODO: non cubic supercells
    n_cells = int(ceil((min_atoms/len(atoms))**(1/3))) or 1

    super_atoms = atoms.repeat((n_cells, n_cells, n_cells))
    super_atoms.info['lattice_constant'] = lattice_constant

    return super_atoms


def a15(lattice_constant, species_A, species_B=None):
    """
    Creates an 8-atom a15-structure in a cubic lattice. A3B.

    Parameters
    ----------
    lattice_constant : float
        Lattice constant for the a15 cubic cell in A.
    species_A : Element
        Atoms to occupy the A sites of the lattice.
    species_B : Element
        Atoms to occupy the B sites of the lattice.
    """

    cell = [[lattice_constant, 0.0, 0.0],
            [0.0, lattice_constant, 0.0],
            [0.0, 0.0, lattice_constant]]

    scaled_positions = [[0.25, 0.50, 0.00],  # A
                        [0.75, 0.50, 0.00],  # A
                        [0.00, 0.25, 0.50],  # A
                        [0.00, 0.75, 0.50],  # A
                        [0.50, 0.00, 0.25],  # A
                        [0.50, 0.00, 0.75],  # A
                        [0.00, 0.00, 0.00],  # B
                        [0.50, 0.50, 0.50]]  # B

    if species_B is None:
        symbols = [species_A]*8
    else:
        symbols = [species_A]*6 + [species_B]*2

    lattice = Atoms(symbols=symbols, scaled_positions=scaled_positions,
                    cell=cell, pbc=True)

    return lattice
