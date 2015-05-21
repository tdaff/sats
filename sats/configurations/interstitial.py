"""
structures.py

Helpers to create different starting environments.

"""

import re

import quippy.structures

from sats.core.elements import Element


def interstitial_dumbbell(lattice, lattice_parameter, species, direction='111',
                          interstitial_species=None, supercell=(1, 1, 1),
                          index=-1):
    """
    Create a crystal structure, such as bcc or fcc, with a dumbbell
    interstitial atom along the given lattice direction.

    Parameters
    ----------
    lattice : str
        Lattice type to generate, valid values are 'bcc' and 'fcc'.
    lattice_parameter : float
        Crystal lattice parameter in A.
    species : int or str
        Identity of species used to construct the BCC lattice.
    direction: str
        Lattice direction along which to create the dumbbell.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    supercell : tuple of int, optional
        Supercell size to place the defect in; if not specified, uses a unit
        cell.
    index : int, optional
        Index of atom with which to create the dumbbell. Added atom is always
        at the end.

    Returns
    -------
    quippy.Atoms
        Lattice with a dumbbell interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    if lattice.lower() == 'bcc':
        bulk = quippy.structures.bcc(a=lattice_parameter,
                                     z=species.atomic_number)
    elif lattice.lower() == 'fcc':
        bulk = quippy.structures.fcc(a=lattice_parameter,
                                     z=species.atomic_number)
    else:
        raise NotImplementedError("Dumbbells not implemented in "
                                  "{0}".format(lattice))

    # catch all False values for supercell
    if supercell:
        bulk = quippy.structures.supercell(bulk, supercell[0], supercell[1],
                                           supercell[2])

    if index < 0:
        index = bulk.n + index

    # parses the digits from any string, ignoring brackets etc.
    direction = [int(x) for x in re.findall('[+-]*[0-9]', direction)]
    displacement = [1.0/x if x else 0.0 for x in direction]
    # Move 1/3 of the <111> atom separation in each direction
    magnitude = lattice_parameter/(12*sum(x**2 for x in displacement))**0.5
    displacement = [x*magnitude for x in displacement]
    bulk.add_atoms(pos=bulk[index].position + displacement,
                   z=interstitial_species.z)
    bulk[index].position -= displacement

    return bulk


def interstitial_crowdion(lattice, lattice_parameter, species,
                          interstitial_species=None, supercell=(1, 1, 1),
                          index=-1):
    """
    Create a crystal structure with a crowdion interstitial atom between
    two close packed atoms.

    Parameters
    ----------
    lattice : str
        Lattice to use for the base crystal structure. Valid choices are 'bcc'
        and 'fcc'.
    lattice_parameter : float
        Lattice parameter in A.
    species : int or str
        Identity of species used to construct the base lattice.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    supercell : tuple of int, optional
        Supercell size to place the defect in; if not specified, uses a unit
        cell.
    index : int, optional
        Index of nearest lattice atom to the crowdion. Additional atom is
        placed along the <111> direction for bcc and <110> direction for fcc.

    Returns
    -------
    quippy.Atoms
        Lattice with a crowdion interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    if lattice.lower() == 'bcc':
        bulk = quippy.structures.bcc(a=lattice_parameter,
                                     z=species.atomic_number)
        displacement = lattice_parameter*0.25
    elif lattice.lower() == 'fcc':
        bulk = quippy.structures.fcc(a=lattice_parameter,
                                     z=species.atomic_number)
        displacement = [lattice_parameter*0.25, lattice_parameter*0.25, 0]
    else:
        raise NotImplementedError("Crowdions not implemented in "
                                  "{0}".format(lattice))

    # catch all False values for supercell
    if supercell:
        bulk = quippy.structures.supercell(bulk, supercell[0], supercell[1],
                                           supercell[2])

    if index < 0:
        index = bulk.n + index

    # Move 0.25 units cells along each axis
    bulk.add_atoms(pos=bulk[index].position + displacement,
                   z=interstitial_species.z)

    return bulk


def interstitial_tetrahedral(lattice, lattice_parameter, species,
                             interstitial_species=None, supercell=(1, 1, 1),
                             index=-1):
    """
    Create a crystal structure with a tetrahedral interstitial atom between
    four atoms.

    Parameters
    ----------
    lattice : str
        Lattice to use for the base crystal structure. Valid choices are 'bcc'
        and 'fcc'.
    lattice_parameter : float
        Lattice parameter in A.
    species : int or str
        Identity of species used to construct the lattice.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    supercell : tuple of int, optional
        Supercell size to place the defect in; if not specified, uses a unit
        cell.
    index : int, optional
        Index of nearest lattice atom to the interstitial. Additional atom is
        placed along the <210> direction in bcc and <221> in fcc.

    Returns
    -------
    quippy.Atoms
        Lattice with a tetrahedral interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    if lattice.lower() == 'bcc':
        bulk = quippy.structures.bcc(a=lattice_parameter,
                                     z=species.atomic_number)
        displacement = [lattice_parameter/4.0, lattice_parameter/2.0, 0]
    elif lattice.lower() == 'fcc':
        bulk = quippy.structures.fcc(a=lattice_parameter,
                                     z=species.atomic_number)
        displacement = [lattice_parameter/4.0,
                        lattice_parameter/4.0,
                        lattice_parameter/2.0]
    else:
        raise NotImplementedError("Tetrahedral defects not implemented in "
                                  "{0}".format(lattice))

    # catch all False values for supercell
    if supercell:
        bulk = quippy.structures.supercell(bulk, supercell[0], supercell[1],
                                           supercell[2])

    if index < 0:
        index = bulk.n + index

    # Move into tetrahedral position
    bulk.add_atoms(pos=bulk[index].position + displacement,
                   z=interstitial_species.z)

    return bulk


def interstitial_octahedral(lattice, lattice_parameter, species,
                            interstitial_species=None, supercell=(1, 1, 1),
                            index=-1):
    """
    Create a crystal structure with an octahedral interstitial atom between
    six atoms.

    Parameters
    ----------
    lattice : str
        Lattice to use for the base crystal structure. Valid choices are 'bcc'
        and 'fcc'.
    lattice_parameter : float
        Lattice parameter in A.
    species : int or str
        Identity of species used to construct the lattice.
    interstitial_species : int or str, optional
        Identity of the interstitial species; if not specified, creates a
        self-interstitial.
    supercell : tuple of int, optional
        Supercell size to place the defect in; if not specified, uses a unit
        cell.
    index : int, optional
        Index of nearest lattice atom to the interstitial. Additional atom is
        placed along the <001> direction.

    Returns
    -------
    quippy.Atoms
        Lattice with an octahedral interstitial atom.
    """

    species = Element(species)
    if interstitial_species is None:
        interstitial_species = species
    else:
        interstitial_species = Element(interstitial_species)

    if lattice.lower() == 'bcc':
        bulk = quippy.structures.bcc(a=lattice_parameter,
                                     z=species.atomic_number)
        displacement = [0, 0, lattice_parameter/2.0]
    elif lattice.lower() == 'fcc':
        bulk = quippy.structures.fcc(a=lattice_parameter,
                                     z=species.atomic_number)
        displacement = [0, 0, lattice_parameter/2.0]
    else:
        raise NotImplementedError("Octahedral defects not implemented in "
                                  "{0}".format(lattice))

    # catch all False values for supercell
    if supercell:
        bulk = quippy.structures.supercell(bulk, supercell[0], supercell[1],
                                           supercell[2])

    if index < 0:
        index = bulk.n + index

    # Move above atom a quarter of a cell
    bulk.add_atoms(pos=bulk[index].position - displacement,
                   z=interstitial_species.z)

    return bulk
