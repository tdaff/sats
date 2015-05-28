"""
interstitial.py

Helpers to create different structures with interstitial defects.

"""

import re

from ase import Atom

from sats.core.elements import Element
from sats.configurations.bulk import bulk


def interstitial(interstitial_id, species, interstitial_species=None,
                 min_atoms=0):
    """
    Helper function to create an interstital based on a string representation
    such as bcc_oct fcc_db_111.

    Parameters
    ----------
    interstitial_id : str
        String representing the interstitial to be created in the form:
        {lattice}_{interstitial}[_{param}].
    species : int or str
        Identity of species used to construct the lattice.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    min_atoms : int, optional
        Construct a supercell to include at least this many atoms.

    Returns
    -------
    interstitial : ase.Atoms
        Lattice with the requested defect.
    """

    irepr = interstitial_id.split('_')
    if 'db' in irepr:
        return interstitial_dumbbell(irepr[0], species, irepr[2],
                                     interstitial_species, min_atoms)
    elif 'crw' in irepr:
        return interstitial_crowdion(irepr[0], species, interstitial_species,
                                     min_atoms)
    elif 'tet' in irepr:
        return interstitial_tetrahedral(irepr[0], species, interstitial_species,
                                        min_atoms)
    elif 'oct':
        return interstitial_octahedral(irepr[0], species, interstitial_species,
                                       min_atoms)
    else:
        raise NotImplementedError("Unknown interstitial {0}.".format(irepr))


def interstitial_dumbbell(lattice, species, direction='111',
                          interstitial_species=None, min_atoms=0, index=-1):
    """
    Create a crystal structure, such as bcc or fcc, with a dumbbell
    interstitial atom along the given lattice direction.

    Parameters
    ----------
    lattice : str
        Lattice type to generate, valid values are 'bcc' and 'fcc' or anything
        accepted by sats bulk.
    species : int or str
        Identity of species used to construct the lattice.
    direction : str
        Lattice direction along which to create the dumbbell.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    min_atoms : int, optional
        Construct a supercell to include at least this many atoms.
    index : int, optional
        Index of atom with which to create the dumbbell. Added atom is always
        at the end.

    Returns
    -------
    ase.Atoms
        Lattice with a dumbbell interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    ibulk = bulk(species, lattice, min_atoms)
    lattice_constant = ibulk.info['lattice_constant']

    if index < 0:
        index += len(ibulk)

    # parses the digits from any string, ignoring brackets etc.
    direction = [int(x) for x in re.findall('[+-]*[0-9]', direction)]
    displacement = [1.0/x if x else 0.0 for x in direction]
    # Move 1/3 of the <111> atom separation in each direction
    magnitude = lattice_constant/(12*sum(x**2 for x in displacement))**0.5
    displacement = [x*magnitude for x in displacement]
    ibulk.append(Atom(interstitial_species,
                      position=ibulk[index].position + displacement))
    ibulk[index].position -= displacement

    return ibulk


def interstitial_crowdion(lattice, species, interstitial_species=None,
                          min_atoms=0, index=-1):
    """
    Create a crystal structure with a crowdion interstitial atom between
    two close packed atoms.

    Parameters
    ----------
    lattice : str
        Lattice to use for the base crystal structure. Valid choices are 'bcc'
        and 'fcc'.
    species : int or str
        Identity of species used to construct the base lattice.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    min_atoms : int, optional
        Construct a supercell to include at least this many atoms.
    index : int, optional
        Index of nearest lattice atom to the crowdion. Additional atom is
        placed along the <111> direction for bcc and <110> direction for fcc.

    Returns
    -------
    ase.Atoms
        Lattice with a crowdion interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    ibulk = bulk(species, lattice, min_atoms)
    lattice_constant = ibulk.info['lattice_constant']

    if lattice.lower() == 'bcc':
        displacement = lattice_constant*0.25
    elif lattice.lower() == 'fcc':
        displacement = [lattice_constant*0.25, lattice_constant*0.25, 0.0]
    elif lattice.lower() == 'hcp':
        displacement = [lattice_constant*0.25, 0.0, 0.0]
    else:
        raise NotImplementedError("Crowdions not implemented in "
                                  "{0}".format(lattice))
    if index < 0:
        index += len(ibulk)

    # Move 0.25 units cells along each axis
    ibulk.append(Atom(interstitial_species,
                      position=ibulk[index].position + displacement))

    return ibulk


def interstitial_tetrahedral(lattice, species, interstitial_species=None,
                             min_atoms=0, index=-1):
    """
    Create a crystal structure with a tetrahedral interstitial atom between
    four atoms.

    Parameters
    ----------
    lattice : str
        Lattice to use for the base crystal structure. Valid choices are 'bcc'
        and 'fcc'.
    species : int or str
        Identity of species used to construct the lattice.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    min_atoms : int, optional
        Construct a supercell to include at least this many atoms.
    index : int, optional
        Index of nearest lattice atom to the interstitial. Additional atom is
        placed along the <210> direction in bcc and <221> in fcc.

    Returns
    -------
    ase.Atoms
        Lattice with a tetrahedral interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    ibulk = bulk(species, lattice, min_atoms)
    lattice_constant = ibulk.info['lattice_constant']

    if lattice.lower() == 'bcc':
        displacement = [lattice_constant/4.0, lattice_constant/2.0, 0]
    elif lattice.lower() == 'fcc':
        displacement = [lattice_constant/4.0,
                        lattice_constant/4.0,
                        lattice_constant/2.0]
    else:
        raise NotImplementedError("Tetrahedral defects not implemented in "
                                  "{0}".format(lattice))

    if index < 0:
        index += len(ibulk)

    # Move into tetrahedral position
    ibulk.append(Atom(interstitial_species,
                      position=ibulk[index].position + displacement))

    return ibulk


def interstitial_octahedral(lattice, species, interstitial_species=None,
                            min_atoms=0, index=-1):
    """
    Create a crystal structure with an octahedral interstitial atom between
    six atoms.

    Parameters
    ----------
    lattice : str
        Lattice to use for the base crystal structure. Valid choices are 'bcc'
        and 'fcc'.
    species : int or str
        Identity of species used to construct the lattice.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    min_atoms : int, optional
        Construct a supercell to include at least this many atoms.
    index : int, optional
        Index of nearest lattice atom to the interstitial. Additional atom is
        placed along the <001> direction.

    Returns
    -------
    ase.Atoms
        Lattice with an octahedral interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    ibulk = bulk(species, lattice, min_atoms)
    lattice_constant = ibulk.info['lattice_constant']

    if lattice.lower() == 'bcc':
        displacement = [0, 0, lattice_constant/2.0]
    elif lattice.lower() == 'fcc':
        displacement = [0, 0, lattice_constant/2.0]
    else:
        raise NotImplementedError("Octahedral defects not implemented in "
                                  "{0}".format(lattice))

    if index < 0:
        index += len(ibulk)

    # Move above atom a quarter of a cell
    ibulk.append(Atom(interstitial_species,
                      position=ibulk[index].position - displacement))

    return ibulk
