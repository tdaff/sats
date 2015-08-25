"""

Bulk material configurations.

"""

from __future__ import division

from math import ceil

from ase import Atoms
from ase.lattice import bulk as ase_bulk

from sats.bulk.rescale import properties
from sats.ui.log import info


def liquid(species, min_atoms=0, supercell=None):
    """
    Generic bcc solid with liquid volume scaling.

    Parameters
    ----------
    species : str or sats.core.elements.Element
        Atomic species used to generate the structures.
    min_atoms : int, optional
        Minimum number of atoms to include in a supercell of the
        bulk. If not specified, you get a unit cell.
    supercell : (int, int, int), optional
        Request a specific supercell of bulk material. If not given,
        max_atoms will be used instead to create a cubic cell.

    Returns
    -------
    liquid : ase.Atoms
        Atoms with liquid equivalent volume (scaled to 0K).

    """

    lattice_constant = properties[species]['lattice_constant']['liquid']

    atoms = ase_bulk(species, 'bcc', a=lattice_constant,
                     orthorhombic=True)

    # TODO: non cubic supercells
    n_cells = int(ceil((min_atoms/len(atoms))**(1/3))) or 1

    if supercell is None:
        info("Generated supercell of {0}x{0}x{0}.".format(n_cells))
        super_atoms = atoms.repeat((n_cells, n_cells, n_cells))
    else:
        info("Generated supercell of {0}.".format(supercell))
        super_atoms = atoms.repeat(supercell)

    info("Structure contains {0} atoms.".format(len(super_atoms)))
    # This ensures that we can find the lattice constant later on
    super_atoms.info['lattice_constant'] = lattice_constant

    return super_atoms


def bulk(species, lattice, min_atoms=0, supercell=None, lattice_constant=None):
    """
    Helper to generate a generic bulk configuration.

    Parameters
    ----------
    species : str or sats.core.elements.Element
        Atomic species used to generate the structures.
    lattice : str
        Bulk lattice type to generate.
    min_atoms : int, optional
        Minimum number of atoms to include in a supercell of the
        bulk. If not specified, you get a unit cell.
    supercell : (int, int, int), optional
        Request a specific supercell of bulk material. If not given,
        max_atoms will be used instead to create a cubic cell.
    lattice_constant : float, optional
        Manually set the lattice constant, rather than using a lookup
        to find the value.

    Returns
    -------
    bulk : ase.Atoms
        Atoms in the desired lattice.
    """

    if lattice_constant is None:
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

    if supercell is None:
        info("Generated supercell of {0}x{0}x{0}.".format(n_cells))
        super_atoms = atoms.repeat((n_cells, n_cells, n_cells))
    else:
        info("Generated supercell of {0}.".format(supercell))
        super_atoms = atoms.repeat(supercell)

    info("Structure contains {0} atoms.".format(len(super_atoms)))
    # This ensures that we can find the lattice constant later on
    super_atoms.info['lattice_constant'] = lattice_constant

    return super_atoms


def primitive(species, lattice, lattice_constant=None):
    """
    Create the primitive unit cell of the given lattice.

    Parameters
    ----------
    species : str or sats.core.elements.Element
        Atomic species used to generate the structure.
    lattice : str
        Lattice type to generate.

    Returns
    -------
    bulk : ase.Atoms
        A primitive cell of the lattice.

    """

    if lattice_constant is None:
        lattice_constant = properties[species]['lattice_constant'][lattice]

    if lattice == 'a15':
        atoms = a15(lattice_constant, species)
    else:
        atoms = ase_bulk(species, lattice, a=lattice_constant,
                         orthorhombic=False)

    # This ensures that we can find the lattice constant later on
    atoms.info['lattice_constant'] = lattice_constant

    return atoms



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
